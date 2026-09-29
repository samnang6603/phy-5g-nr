#include <cmath>
#include <complex>
#include <cstddef>
#include "generator.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../cluster_profiles/predefined_profiles.hpp"
#include "../ray_coupling/coupling.hpp"
#include "../antenna_array/antenna.hpp"
#include "../mobility/mobility.hpp"

static constexpr float FULL_CIRCLE_DEGREES = 360.0f;
static constexpr float HALF_CIRCLE_DEGREES = 180.0f;
static constexpr std::size_t NUM_FIELD_SPHERICAL_ANGLES = 2; // 2 for F theta and F phi

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
    struct Dim {
        std::size_t T = 1;
        std::size_t X = NULL;
        std::size_t M;
        std::size_t P;
        std::size_t R;
    };

    /***************** Function Declarations ********************/
    static std::vector<float> compute_LoS_cluster_gain(
        Dim& D,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const std::vector<float>& phi,
        const std::vector<std::size_t>& coupling,
        const float XPR
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

        Dim D{1,
            0,
            NUMBER_OF_RAYS_PER_CLUSTER,
            numTx,
            numRx,                    
            };

        const std::vector<PropagationCondition>& cluster_types = pdp_conf.ClusterTypes;
        const std::size_t L = cluster_types.size();
        const float XPR = pdp_conf.XPR;

        if (hasLoS) {
            std::vector<float> los_gain = compute_LoS_cluster_gain(
                D, 
                ant_sys_conf, 
                pdp_conf, 
                phi, 
                coupling, 
                XPR
            );
        }



    }

    static std::vector<float> compute_LoS_cluster_gain(
        Dim& D,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const std::vector<float>& phi,
        const std::vector<std::size_t>& coupling,
        const float XPR
    ) {
        // Compute cluster gain for LoS component

        D.X = 0;

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

        // ----------------------------------Processing for Tx-----------------------------------------------------------------
        std::array<float, antenna::geometry::NUM_3D_AXIS> rhat_tx = 
            antenna::geometry::get_LoS_spherical_unit_vector(phi_AoD, theta_ZoD);

        // Get transmit antenna/subarray location vector
        const antenna::AntennaArrayConfig& tx_ant_conf = ant_sys_conf.TransmitAntennaArray;
        const std::vector<float>& tx_radiator_pos = tx_ant_conf.FieldPattern.ElementPositions;
        const antenna::geometry::PositionConfig& pos_conf = tx_ant_conf.Position;

        // Calculate the location vector dbar of Tx
        std::vector<float> dbar_tx(tx_radiator_pos.size());
        for (std::size_t i = 0; i < tx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
            dbar_tx[i] = tx_radiator_pos[i]*lambda_0 + pos_conf.x;
            dbar_tx[i + 1] = tx_radiator_pos[i + 1]*lambda_0 + pos_conf.y;
            dbar_tx[i + 2] = tx_radiator_pos[i + 2]*lambda_0 + pos_conf.z;
        }
        
        // Allocation for Tx field term and location term
        const std::size_t numTx = ant_sys_conf.NumInputSignals;
        const std::size_t term_size = NUMBER_OF_RAYS_PER_CLUSTER*numTx;
        std::vector<std::complex<float>> tx_field_term(NUM_FIELD_SPHERICAL_ANGLES*term_size); 
        std::vector<std::complex<float>> tx_loc_term(term_size);

        //std::vector<float> field_term_tx(numTx*antenna::NUM_MAX_POLARIZATION);

        // Compute 
        for (std::size_t s = 0; s < numTx; ++s) {

            auto field_term_tx = antenna::field_pattern::compute_LoS_field_term(
                tx_ant_conf,
                theta_ZoD, 
                phi_AoD, 
                s
            );

            auto loc_term_tx = antenna::field_pattern::get_LoS_location_term(
                rhat_tx,
                dbar_tx,
                lambda_0,
                s
            );

            //std::cout << "Aha!" << std::endl;
        }

        //-----------------------------------End Processing for Tx-----------------------------------------------------------------

        const std::size_t numRx = ant_sys_conf.NumOutputSignals;





        std::vector<float> x(2,0.0f);
        return x;



        /*
        // Get corresponding ray coupling
        const std::size_t* coup_ptr = coupling.data();

        // Get corresponding initial phase
        const float* phi_ptr = phi.data();

        // Per cluster parameter
        const float C_ASD = pdp_conf.AngleSpreads.C_ASD;
        const float C_ASA = pdp_conf.AngleSpreads.C_ASA;
        const float C_ZSD = pdp_conf.AngleSpreads.C_ZSD;
        const float C_ZSA = pdp_conf.AngleSpreads.C_ZSA; */

        
        

    }

}