#include <cmath>
#include <complex>
#include <cstddef>
#include "generator.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../cluster_profiles/predefined_profiles.hpp"
#include "../ray_coupling/coupling.hpp"
#include "../antenna_array/antenna.hpp"
#include "../mobility/mobility.hpp"

/************************* Constants *****************************/
static constexpr float FULL_CIRCLE_DEGREES = 360.0f;
static constexpr float HALF_CIRCLE_DEGREES = 180.0f;
static constexpr std::size_t NUM_FIELD_SPHERICAL_ANGLES = 2; // 2 for F theta and F phi
static constexpr std::size_t NUM_PATH_ANGLES = 4;

// MATLAB-compatible mod(x, y) for positive y.
static inline constexpr float MMOD(float x, float y) {
    return std::fmod(std::fmod(x,y) + y, y);
}

// Wrap azimuth angles to [-180, 180]
// Shift by 180, wrap modulo 360, then center at 0
static inline constexpr float WRAP_AZIMUTH_ANGLE_180(float angle) {
    return MMOD(angle + HALF_CIRCLE_DEGREES, FULL_CIRCLE_DEGREES) - HALF_CIRCLE_DEGREES;
}

// Wrap zenith angle to [0, 360), then reflect (180, 360) into (180, 0).
// Result is in [0, 180].
static inline constexpr float WRAP_ZENITH_ANGLE(float angle) {
    float angle2 = MMOD(angle, FULL_CIRCLE_DEGREES);
    return angle2 > HALF_CIRCLE_DEGREES ? FULL_CIRCLE_DEGREES - angle2 : angle2;
}


namespace channels::cdl::response {

    /***************** Type Definitions ************************/
    struct DimConfig {
        std::size_t T = 1;
        std::size_t X = 1;
        std::size_t M = 1;
        std::size_t P = 1;
        std::size_t R = 1;

        std::size_t NumDim() const noexcept {
            return T*X*M*P*R;
        }
    };

    /***************** Function Declarations ********************/
    static std::vector<std::complex<float>> compute_LoS_cluster_gain(
        const DimConfig& dim_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const std::vector<float>& phi,
        const std::vector<std::size_t>& coupling,
        const float XPR
    );

    static std::vector<std::complex<float>> calculate_polarization_matrix(
        const std::vector<float>& phi,
        const float XPR,
        const std::vector<PropagationCondition> cluster_type,
        const PropagationCondition this_cluster_type
    );
    
    /***************** Function Implementations ********************/
    std::vector<float> generate_static_path_gains(
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& phi
    ) {
        // Generate static cdl path gains without Doppler effect

        const std::size_t numTx = ant_sys_conf.NumInputSignals;
        const std::size_t numRx = ant_sys_conf.NumOutputSignals;
        const bool hasLoS = pdp_conf.HasLoS;

        DimConfig D{1,
                    1,
                    NUMBER_OF_RAYS_PER_CLUSTER,
                    numTx,
                    numRx,                    
                };

        const std::vector<PropagationCondition>& cluster_types = pdp_conf.ClusterTypes;
        const std::size_t L = cluster_types.size();
        const float XPR = pdp_conf.XPR;

        if (hasLoS) {
            std::vector<std::complex<float>> Hstatic_los = compute_LoS_cluster_gain(
                D, 
                ant_sys_conf, 
                pdp_conf, 
                phi, 
                coupling, 
                XPR
            );
        }

        std::vector<float> x(5,1.0f);

        return x;

    }

    static std::vector<std::complex<float>> compute_LoS_cluster_gain(
        const DimConfig& dim_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const std::vector<float>& phi,
        const std::vector<std::size_t>& coupling,
        const float XPR
    ) {
        // Compute cluster gain for LoS component

        const pdp::CDLCluster pdp_LoS_cluster = pdp_conf.Table[0];
        const float pow_linear = std::powf(10.0f, pdp_LoS_cluster.power_db/10.0f);
        
        // 7.5 Equation 7.5-13 for AoA and AoD (ray_offset = 0 for LoS)
        // Also wrap AoA and AoD to [-180, 180]
        const float phi_AoA = WRAP_AZIMUTH_ANGLE_180(pdp_LoS_cluster.aoa_deg);
        const float phi_AoD = WRAP_AZIMUTH_ANGLE_180(pdp_LoS_cluster.aod_deg);
        
        // 7.5 Equation 7.5-18 for ZoA and ZoD (ray_offset = 0 for LoS)
        const float theta_ZoA = WRAP_ZENITH_ANGLE(pdp_LoS_cluster.zoa_deg);
        const float theta_ZoD = WRAP_ZENITH_ANGLE(pdp_LoS_cluster.zod_deg);
        
        // TODO: Angle scaling section 7.7.5.1 to be implemented

        const float lambda_0 = ant_sys_conf.Wavelength();
        const float lambda_0_inv = 1/lambda_0;

        // -------------Processing for Tx----------------------------------------------------------------------------------------
        std::array<float, antenna::geometry::NUM_3D_AXIS> rhat_tx = 
            antenna::geometry::get_LoS_spherical_unit_vector(phi_AoD, theta_ZoD); // angles of departure

        // Get Tx antenna/subarray location vector
        const antenna::AntennaArrayConfig& tx_ant_conf = ant_sys_conf.TransmitAntennaArray;
        const std::vector<float>& tx_radiator_pos = tx_ant_conf.FieldPattern.ElementPositions;
        const antenna::geometry::PositionConfig& tx_pos_conf = tx_ant_conf.Position;

        // Calculate the location vector dbar of Tx
        std::vector<float> dbar_tx(tx_radiator_pos.size());
        for (std::size_t i = 0; i < tx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
            dbar_tx[i] = tx_radiator_pos[i]*lambda_0 + tx_pos_conf.x;
            dbar_tx[i + 1] = tx_radiator_pos[i + 1]*lambda_0 + tx_pos_conf.y;
            dbar_tx[i + 2] = tx_radiator_pos[i + 2]*lambda_0 + tx_pos_conf.z;
        }
        
        // Allocation for Tx field term and location term
        const std::size_t numTx = ant_sys_conf.NumInputSignals;
        std::vector<float> tx_field_term(numTx*NUM_FIELD_SPHERICAL_ANGLES);
        std::vector<std::complex<float>> tx_loc_term(numTx*NUMBER_OF_RAYS_PER_CLUSTER);

        // Compute field and location terms
        for (std::size_t s = 0; s < numTx; ++s) {

            float* field_term = tx_field_term.data() + s*NUM_FIELD_SPHERICAL_ANGLES;

            antenna::field_pattern::compute_LoS_field_term(
                field_term,
                tx_ant_conf,
                theta_ZoD, 
                phi_AoD, 
                s
            );

            std::complex<float>* loc_term = tx_loc_term.data() + s*NUMBER_OF_RAYS_PER_CLUSTER;

            antenna::field_pattern::get_LoS_location_term(
                loc_term,
                rhat_tx,
                dbar_tx,
                lambda_0_inv,
                s
            );
        }
        //------------End Processing for Tx----------------------------------------------------------------------------------------


        //--------------Processing for Rx------------------------------------------------------------------------------------------
        std::array<float, antenna::geometry::NUM_3D_AXIS> rhat_rx = 
            antenna::geometry::get_LoS_spherical_unit_vector(phi_AoA, theta_ZoA); // angles of arrival

        // Get Rx antenna/subarray location vector
        const antenna::AntennaArrayConfig&rx_ant_conf = ant_sys_conf.ReceiveAntennaArray;
        const std::vector<float>& rx_radiator_pos = rx_ant_conf.FieldPattern.ElementPositions;
        const antenna::geometry::PositionConfig& rx_pos_conf = rx_ant_conf.Position;

        // Calculate the location vector dbar of Tx
        std::vector<float> dbar_rx(rx_radiator_pos.size());
        for (std::size_t i = 0; i < rx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
            dbar_rx[i] = rx_radiator_pos[i]*lambda_0 + rx_pos_conf.x;
            dbar_rx[i + 1] = rx_radiator_pos[i + 1]*lambda_0 + rx_pos_conf.y;
            dbar_rx[i + 2] = rx_radiator_pos[i + 2]*lambda_0 + rx_pos_conf.z;
        }

        // Allocation for Tx field term and location term
        const std::size_t numRx = ant_sys_conf.NumOutputSignals;
        std::vector<float> rx_field_term(numRx*NUM_FIELD_SPHERICAL_ANGLES);
        std::vector<std::complex<float>> rx_loc_term(numRx*NUMBER_OF_RAYS_PER_CLUSTER);

        // Compute field and location terms
        for (std::size_t u = 0; u < numRx; ++u) {

            float* field_term = rx_field_term.data() + u*NUM_FIELD_SPHERICAL_ANGLES;

            antenna::field_pattern::compute_LoS_field_term(
                field_term,
                rx_ant_conf,
                theta_ZoA, 
                phi_AoA, 
                u
            );

            std::complex<float>* loc_term = rx_loc_term.data() + u*NUMBER_OF_RAYS_PER_CLUSTER;

            antenna::field_pattern::get_LoS_location_term(
                loc_term,
                rhat_rx,
                dbar_rx,
                lambda_0_inv,
                u
            );
        }
        //------------End Processing for Rx----------------------------------------------------------------------------------------
        auto polterm = calculate_polarization_matrix(
            phi, 
            XPR, 
            pdp_conf.ClusterTypes, 
            PropagationCondition::LOS
        );

        
        // Calculate the MONSTROUS Equation 7.5-28 excluding Doppler (the
        // last term). Doppler is to be calculated in time-varying channel response
        // First, compute all the field terms combined. In 7.5-28, all the 
        // field terms are combined as below:
        /*
            (thth): theta-theta, (thph): theta-phi
            ksqrinv: 1/sqrt(kappa)

        fieldTerms =

                 T
        [F_rx_th]   [exp(jPhi_thth)          ksqrinv*exp(jPhi_thph]   [F_tx_th]
        |       | * |                                             | * |       |
        [F_rx_ph]   [ksqrinv*exp(jPhi_phth)         exp(jPhi_phph)]   [F_rx_ph]
             
        */
        std::vector<std::complex<float>> Hstatic(numRx*numTx); // no rays consideration because this is LoS
        std::complex<float>* h = Hstatic.data();
        // This implmentation doesn't include matrix multiplication as in 7.5-28 but an expanded
        // vectorized version for optimization.
        // For this LoS, we really don't need to iterate across the rays, since all rays have uniform equal power.
        // All it matters are the polarization matrix (polterm)
        for (std::size_t u = 0; u < numRx; ++u) {

            for (std::size_t s = 0; s < numTx; ++s) {

                const std::size_t u_stride = u*NUM_FIELD_SPHERICAL_ANGLES;
                const std::size_t s_stride = s*NUM_FIELD_SPHERICAL_ANGLES;

                const float* rxft1 = rx_field_term.data() + u_stride;
                const float* rxft2 = rx_field_term.data() + u_stride + 1;

                const float* txft1 = tx_field_term.data() + s_stride;
                const float* txft2 = tx_field_term.data() + s_stride + 1;

                const std::complex<float> tmp0 = (*rxft1)*polterm[0];
                const std::complex<float> tmp1 = (*rxft2)*polterm[1];
                const std::complex<float> tmp2 = (*rxft1)*polterm[2];
                const std::complex<float> tmp3 = (*rxft2)*polterm[3];
                const std::complex<float> sum_field_terms = (tmp0 + tmp1)*(*txft1) + (tmp2 + tmp3)*(*txft2);

                /*
                  Then, compute all the location terms combined. In 7.5-28, all the 
                  location terms are combined as below (excluding Doppler term)
                
                                    T                            T
                             /     r_rx*dbar_rx  \        /     r_tx*dbar_tx  \    
                         exp| j2pi--------------- | * exp| j2pi--------------- |
                             \       lambda_0    /        \       lambda_0    / 
            
                */

                const std::complex<float>* rxloc = rx_loc_term.data() + u*NUMBER_OF_RAYS_PER_CLUSTER;
                const std::complex<float>* txloc = tx_loc_term.data() + s*NUMBER_OF_RAYS_PER_CLUSTER;

                *h++ = (*rxloc++)*(*txloc++)*sum_field_terms;
            }
        }
        return Hstatic;
    }


    static std::vector<std::complex<float>> calculate_polarization_matrix(
        const std::vector<float>& phi,
        const float XPR,
        const std::vector<PropagationCondition> cluster_type,
        const PropagationCondition this_cluster_type
    ) {
        // Calculate polarization matrix for the specified cluster type
        // The Phi (Angular interaction) is arranged like this:
        // | theta-theta |
        // |   phi-theta |
        // | theta-phi   |
        // |   phi-phi   |

        std::size_t num_cluster = cluster_type.size();
        std::vector<std::complex<float>> cphi(num_cluster);
        std::size_t count = 0;
        for (std::size_t i = 0; i < num_cluster; ++i) {
            if (cluster_type[i] == this_cluster_type) {
                cphi[i] = std::complex{
                    std::cos(phi[i]), 
                    std::sin(phi[i])
                };
                ++count;
            }
        }

        std::vector<std::complex<float>> polmat(count*NUM_PATH_ANGLES);

        switch (this_cluster_type) {
            case PropagationCondition::LOS:

                // For LoS case, it is really simple
                // only theta-theta and phi-phi matter, no cross angle
                polmat[0] = cphi[0];
                polmat[3] = -cphi[0]; // phi-phi is negative

            case PropagationCondition::NLOS:


            default:

        }

        return polmat;

    }

} // namespace channels::cdl::response