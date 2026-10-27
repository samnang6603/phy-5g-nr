#include <cmath>
#include <complex>
#include <vector>
#include <span>
#include <cstddef>
#include "generator.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../cluster_profiles/predefined_profiles.hpp"
#include "../antenna_array/antenna.hpp"

/************************* Constants *****************************/
static constexpr float FULL_CIRCLE_DEGREES = 360.0f;
static constexpr float HALF_CIRCLE_DEGREES = 180.0f;
static constexpr std::size_t NUM_FIELD_SPHERICAL_ANGLES = 2; // 2 for F theta and F phi
static constexpr std::size_t NUM_PATH_ANGLES = 4;

// Wrap azimuth angles to [-180, 180]
// Shift by 180, wrap modulo 360, then center at 0
static inline constexpr float WRAP_AZIMUTH_ANGLE_180(const float angle)
{
    return angle - FULL_CIRCLE_DEGREES*std::floor((angle + HALF_CIRCLE_DEGREES)/FULL_CIRCLE_DEGREES);
}

// Wrap zenith angle to [0, 360), then reflect (180, 360) into (180, 0).
// Result is in [0, 180].
static inline float WRAP_ZENITH_ANGLE(const float angle)
{
    const float x = angle - FULL_CIRCLE_DEGREES*std::floor(angle/FULL_CIRCLE_DEGREES);
    return HALF_CIRCLE_DEGREES - std::abs(x - HALF_CIRCLE_DEGREES);
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
    namespace nlos {

        static void compute_cluster_gain(
            StaticPathGainConfig& nlos_gain_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const pdp::DelayProfileConfig& pdp_conf,
            const std::vector<float>& Phi,
            const std::vector<std::size_t>& coupling
        );
    } // namespace nlos

    namespace los {

        static void compute_cluster_gain(
            StaticPathGainConfig& los_gain_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const pdp::DelayProfileConfig& pdp_conf,
            const std::vector<float>& Phi
        );
    } //namespace los

    static std::vector<std::complex<float>> calculate_polarization_matrix(
        const std::vector<float>& Phi,
        const float XPR,
        const std::vector<PropagationCondition> all_cluster_type,
        const PropagationCondition this_cluster_type
    );

    void generate_static_path_gains(
        StaticPathGainConfig& los_gain_conf,
        StaticPathGainConfig& nlos_gain_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& Phi
    ) {
        // Generate static cdl path gains without Doppler effect

        const std::size_t num_Tx = ant_sys_conf.NumInputSignals;
        const std::size_t num_Rx = ant_sys_conf.NumOutputSignals;
        const bool hasLoS = pdp_conf.HasLoS;

        DimConfig D{1,
                    1,
                    NUM_RAYS_PER_CLUSTER,
                    num_Tx,
                    num_Rx,                    
                };

        if (hasLoS) {
            los::compute_cluster_gain(
                los_gain_conf,
                channel_ctrl_conf,
                ant_sys_conf,
                pdp_conf,
                Phi
            );
        }

        nlos::compute_cluster_gain(
            nlos_gain_conf,
            channel_ctrl_conf,
            ant_sys_conf, 
            pdp_conf, 
            Phi, 
            coupling
        );

        // TODO: Subcluster gain

    }

    namespace nlos {

        static void compute_cluster_gain(
            StaticPathGainConfig& nlos_gain_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const pdp::DelayProfileConfig& pdp_conf,
            const std::vector<float>& Phi,
            const std::vector<std::size_t>& coupling
        ) {
            // Compute cluster gain for NLoS components.

            const std::size_t los_ind_start = static_cast<std::size_t>(pdp_conf.HasLoS);
            const std::size_t nlos_ind_start = los_ind_start ? 1 : 0;
            const std::size_t num_all_cluster = pdp_conf.ClusterTypes.size();
            const std::size_t num_nlos_cluster =  num_all_cluster - static_cast<std::size_t>(pdp_conf.HasLoS);

            const std::span<const pdp::CDLCluster> pdp_NLoS_cluster =
                std::span(pdp_conf.Table).subspan(nlos_ind_start);

            const pdp::AngleSpreadsList& angle_spread = pdp_conf.AngleSpreads;
            const float C_ASD = angle_spread.C_ASD;
            const float C_ASA = angle_spread.C_ASA;
            const float C_ZSD = angle_spread.C_ZSD;
            const float C_ZSA = angle_spread.C_ZSA;

            // Calculate angle parameters
            std::size_t total_cluster_samples = num_all_cluster*NUM_RAYS_PER_CLUSTER;
            std::size_t total_nlos_samples = num_nlos_cluster*NUM_RAYS_PER_CLUSTER;
            std::vector<float> phi_AoA(total_nlos_samples);
            float* phi_aoa_ptr = phi_AoA.data();

            std::vector<float> phi_AoD(total_nlos_samples);
            float* phi_aod_ptr = phi_AoD.data();
            
            std::vector<float> theta_ZoA(total_nlos_samples);
            float* theta_zoa_ptr = theta_ZoA.data();

            std::vector<float> theta_ZoD(total_nlos_samples);
            float* theta_zod_ptr = theta_ZoD.data();

            for (const auto& r : RAY_OFFSET_ALPHA) {

                for (const auto& v : pdp_NLoS_cluster) {
                    // 7.5 Eq. 7.5-13 for AoA and AoD
                    // Also wrap AoA and AoD to [-180, 180]
                    *phi_aoa_ptr++ = WRAP_AZIMUTH_ANGLE_180(v.aoa_deg + C_ASA*r);
                    *phi_aod_ptr++ = WRAP_AZIMUTH_ANGLE_180(v.aod_deg + C_ASD*r);

                    // 7.5 Eq. 7.5-18 for ZoA and ZoD
                    *theta_zoa_ptr++ = WRAP_ZENITH_ANGLE(v.zoa_deg + C_ZSA*r);
                    *theta_zod_ptr++ = WRAP_ZENITH_ANGLE(v.zod_deg + C_ZSD*r);
                }
            }

            // Ray coupling angle parameters indices swap
            // Allocate the angles
            std::vector<float> phi_AoA_rc(total_nlos_samples);
            float* phi_aoa_rc_ptr = phi_AoA_rc.data();
            
            std::vector<float> theta_ZoA_rc(total_nlos_samples);
            float* theta_zoa_rc_ptr = theta_ZoA_rc.data();

            std::vector<float> theta_ZoD_rc(total_nlos_samples);
            float* theta_zod_rc_ptr = theta_ZoD_rc.data();

            std::size_t aod_stride = total_cluster_samples*2;
            std::size_t zod_stride = total_cluster_samples;

            std::size_t coupling_idx;
            for (std::size_t r = 0; r < NUM_RAYS_PER_CLUSTER; ++r) {

                for (std::size_t n = 0; n < num_nlos_cluster; ++n) {

                    // Coupling permutations
                    // 1. AoD -> AoA coupling
                    // 2: AoD -> ZoA coupling
                    // 3. AoD -> ZoD coupling
                    const std::size_t ind_start = los_ind_start + r + n + num_nlos_cluster*r;
                    //std::size_t ind_start = n + num_nlos_cluster*r;
                    
                    // AoD -> AoA coupling
                    coupling_idx = coupling[ind_start];
                    *phi_aoa_rc_ptr++ = phi_AoA[coupling_idx];
                    
                    // 7.5 Eq. 7.5-18 for ZoA and ZoD
                    // AoD -> ZoD coupling
                    coupling_idx = coupling[ind_start + aod_stride];
                    *theta_zod_rc_ptr++ = theta_ZoD[coupling_idx];

                    // TODO: to be removed when ZoD -> ZoA coupling + rearrange in AoD order combining done
                    // ZoD -> ZoA coupling
                    coupling_idx = coupling[ind_start + zod_stride];
                    *theta_zoa_rc_ptr++ = theta_ZoA[coupling_idx]; 

                    /* TODO: combine ZoD -> ZoA coupling + rearrange in AoD order
                    // ZoD -> ZoA coupling, then rearrange in AoD order
                    const std::size_t zoa_idx = coupling[ind_start + aod_stride];

                    // Convert the intermediate index to the flat coupling array
                    const std::size_t zoa_coupling_idx = coupling[los_ind_start + zoa_idx + zod_stride];

                    *theta_zoa_rc_ptr++ = theta_ZoA[zoa_coupling_idx]; 
                    */
                }
            }

            // TODO: to be removed when ZoD -> ZoA coupling + rearrange in AoD order combining done
            // Rearrange ZoA in AoD order
            std::vector<float> theta_ZoA_rc2(total_nlos_samples);
            theta_zoa_rc_ptr = theta_ZoA_rc2.data();
            for (std::size_t r = 0; r < NUM_RAYS_PER_CLUSTER; ++r) {

                for (std::size_t n = 0; n < num_nlos_cluster; ++n) {

                    const std::size_t ind_start = los_ind_start + r + n + num_nlos_cluster*r;
                    coupling_idx = coupling[ind_start + aod_stride];
                    *theta_zoa_rc_ptr++ = theta_ZoA_rc[coupling_idx];

                }
            } 

            // TODO: Angle scaling section 7.7.5.1 to be implemented

            // Get wavelength
            const float lambda_0 = ant_sys_conf.Wavelength();
            const float lambda_0_inv = 1/lambda_0;

            // -------------Processing for Tx----------------------------------------------------------------------------------------
            const antenna::geometry::nlos::SphericalDirections tx_sph_dir_conf = 
                antenna::geometry::nlos::get_spherical_unit_vectors(phi_AoD, theta_ZoD_rc); // angles of departure
            
            nlos_gain_conf.TxSphericalUnitVector.assign(
                tx_sph_dir_conf.rhat.begin(),
                tx_sph_dir_conf.rhat.end()
            );

            // Get Tx antenna/subarray location vector
            const antenna::AntennaArrayConfig& tx_ant_conf = ant_sys_conf.TransmitAntennaArray;
            const std::vector<float>& tx_radiator_pos = tx_ant_conf.FieldPattern.ElementPositions;
            const antenna::geometry::PositionConfig& tx_pos_conf = tx_ant_conf.Position;

            // Calculate the location vector dbar of Tx
            std::vector<float> tx_dbar(tx_radiator_pos.size());
            for (std::size_t i = 0; i < tx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
                tx_dbar[i] = tx_radiator_pos[i]*lambda_0 + tx_pos_conf.x;
                tx_dbar[i + 1] = tx_radiator_pos[i + 1]*lambda_0 + tx_pos_conf.y;
                tx_dbar[i + 2] = tx_radiator_pos[i + 2]*lambda_0 + tx_pos_conf.z;
            }

            // Allocation for Tx field term and location term
            const std::size_t num_Tx = ant_sys_conf.NumInputSignals;
            const std::size_t tx_term_size = num_Tx*num_nlos_cluster*NUM_RAYS_PER_CLUSTER;
            std::vector<float> tx_field_term(tx_term_size*NUM_FIELD_SPHERICAL_ANGLES);
            std::vector<std::complex<float>> tx_loc_term(tx_term_size);

            for (std::size_t s = 0; s < num_Tx; ++s) {

                float* field_term = tx_field_term.data() + s*NUM_FIELD_SPHERICAL_ANGLES*total_nlos_samples;

                antenna::field_pattern::nlos::compute_field_term(
                    field_term,
                    tx_ant_conf,
                    tx_sph_dir_conf,
                    s
                ); 

                std::complex<float>* loc_term = tx_loc_term.data() + s*total_nlos_samples;

                antenna::field_pattern::nlos::get_location_term(
                    loc_term,
                    tx_sph_dir_conf.rhat,
                    tx_dbar,
                    total_nlos_samples,
                    lambda_0_inv,
                    s
                );
            }
            //------------End Processing for Tx----------------------------------------------------------------------------------------

            //--------------Processing for Rx------------------------------------------------------------------------------------------
            const antenna::geometry::nlos::SphericalDirections rx_sph_dir_conf = 
                antenna::geometry::nlos::get_spherical_unit_vectors(phi_AoA_rc, theta_ZoA_rc2); // angles of arrival
            
            nlos_gain_conf.RxSphericalUnitVector.assign(
                rx_sph_dir_conf.rhat.begin(),
                rx_sph_dir_conf.rhat.end()
            );

            // Get Rx antenna/subarray location vector
            const antenna::AntennaArrayConfig& rx_ant_conf = ant_sys_conf.ReceiveAntennaArray;
            const std::vector<float>& rx_radiator_pos = rx_ant_conf.FieldPattern.ElementPositions;
            const antenna::geometry::PositionConfig& rx_pos_conf = rx_ant_conf.Position;

            // Calculate the location vector dbar of Rx
            std::vector<float> rx_dbar(rx_radiator_pos.size());
            for (std::size_t i = 0; i < rx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
                rx_dbar[i] = rx_radiator_pos[i]*lambda_0 + rx_pos_conf.x;
                rx_dbar[i + 1] = rx_radiator_pos[i + 1]*lambda_0 + rx_pos_conf.y;
                rx_dbar[i + 2] = rx_radiator_pos[i + 2]*lambda_0 + rx_pos_conf.z;
            }

            // Allocation for Tx field term and location term
            const std::size_t num_Rx = ant_sys_conf.NumOutputSignals;
            const std::size_t rx_term_size = num_Rx*num_nlos_cluster*NUM_RAYS_PER_CLUSTER;
            std::vector<float> rx_field_term(rx_term_size*NUM_FIELD_SPHERICAL_ANGLES);
            std::vector<std::complex<float>> rx_loc_term(rx_term_size);

            for (std::size_t u = 0; u < num_Rx; ++u) {

                float* field_term = rx_field_term.data() + u*NUM_FIELD_SPHERICAL_ANGLES*total_nlos_samples;

                antenna::field_pattern::nlos::compute_field_term(
                    field_term,
                    rx_ant_conf,
                    rx_sph_dir_conf,
                    u
                ); 

                std::complex<float>* loc_term = rx_loc_term.data() + u*total_nlos_samples;

                antenna::field_pattern::nlos::get_location_term(
                    loc_term,
                    rx_sph_dir_conf.rhat,
                    rx_dbar,
                    total_nlos_samples,
                    lambda_0_inv,
                    u
                );
            }
            //------------End Processing for Rx----------------------------------------------------------------------------------------

            // Calculate polarization matrix
            const float XPR = pdp_conf.XPR;
            std::vector<std::complex<float>> polterm = calculate_polarization_matrix(
                Phi, 
                XPR, 
                pdp_conf.ClusterTypes, 
                PropagationCondition::NLOS
            );

            // Calculate output normalization gain
            const float norm_output = channel_ctrl_conf.NormalizeChannelOutput ? 1.0f/std::sqrt(num_Rx) : 1.0f;

            // Get active ray
            // TODO: find out if finding active ray is truly needed once more features are added like subclustering
             
            // Scaling factor: ksqrtinv in 7.5-28
            // Expand ksqrtinv from num_nlos_cluster to total_nlos_samples
            // (num_nlos_cluster * NUM_RAYS_PER_CLUSTER) to avoid index bookkeeping
            // and enable contiguous scaling-factor reads when computing H_static.
            std::vector<float> ksqrtinv(total_nlos_samples);
            const bool normalize = channel_ctrl_conf.NormalizePathGains;
            float avg_path_gains = 0.0f;
            for (std::size_t i = 0; i < num_nlos_cluster; ++i) {

                // Compute amplitude once per cluster.
                const float amplitude = std::pow(10.0f, pdp_NLoS_cluster[i].power_db*0.05f);

                // Accumulate linear cluster power (before ray normalization).
                if (normalize) {
                    avg_path_gains += amplitude*amplitude;
                }

                // Compute the scaling factor once per cluster and reuse it for all rays.
                const float scale = amplitude*INV_SQRT_NUM_RAYS_PER_CLUSTER;

                // Store in cluster-fastest order, matching H_static's layout
                // Writes are strided by num_nlos_cluster, but this small initialization
                // cost is preferable to recomputing std::pow for every ray.
                for (std::size_t r = 0; r < NUM_RAYS_PER_CLUSTER; ++r) {
                    ksqrtinv[r * num_nlos_cluster + i] = scale;
                }
            }
            float norm_path_gain = 1.0f;
            if (normalize) {

                if (pdp_conf.HasLoS) {

                    avg_path_gains += std::pow(10.0f, pdp_conf.Table[0].power_db*0.1f);
                }

                norm_path_gain = 1.0f/std::sqrt(avg_path_gains);
            }

            const float norm_path_and_output = norm_path_gain*norm_output;

            // Calculate the MONSTROUS Eq. 7.5-28 excluding Doppler (the
            // last term). Doppler is to be calculated in time-varying channel response
            nlos_gain_conf.RayResponse.resize(total_nlos_samples*num_Rx*num_Tx);
            std::complex<float>* h = nlos_gain_conf.RayResponse.data();
            std::complex<float>* polterm_ptr;
            const std::size_t field_stride = NUM_FIELD_SPHERICAL_ANGLES*total_nlos_samples;
            
            for (std::size_t u = 0; u < num_Rx; ++u) {

                const float* rxft_base = rx_field_term.data() + u*field_stride;
                const std::complex<float>* rxloc_base = rx_loc_term.data() + u * total_nlos_samples;

                for (std::size_t s = 0; s < num_Tx; ++s) {

                    // Reset to the beginning of this receiver's data
                    const float* rxft = rxft_base;
                    const std::complex<float>* rxloc = rxloc_base;

                    const float* txft = tx_field_term.data() + s*field_stride;
                    const std::complex<float>* txloc = tx_loc_term.data() + s*total_nlos_samples;

                    polterm_ptr = polterm.data();

                    for (std::size_t i = 0; i < total_nlos_samples; ++i) {
                        // First, compute all the field terms combined. In 7.5-28, all the 
                        // field terms are combined as below:
                        /*
                            (thth): theta-theta, (thph): theta-phi
                            ksqrtinv: 1/sqrt(kappa)

                         fieldTerms =

                                  T
                         [F_rx_th]   [exp(jPhi_thth)          ksqrinv*exp(jPhi_thph]   [F_tx_th]
                         |       | * |                                             | * |       |
                         [F_rx_ph]   [ksqrinv*exp(jPhi_phth)         exp(jPhi_phph)]   [F_rx_ph]
                       
                        */

                        // This implmentation doesn't include matrix multiplication as in 7.5-28 
                        // but an expanded vectorized version for optimization.
                        // ksqrtinv included later
                        const std::complex<float> tmp0 = (rxft[0])*(*polterm_ptr++);  // [0] is F_theta
                        const std::complex<float> tmp1 = (rxft[1])*(*polterm_ptr++);  // [1] is F_phi
                        const std::complex<float> tmp2 = (rxft[0])*(*polterm_ptr++);
                        const std::complex<float> tmp3 = (rxft[1])*(*polterm_ptr++);
                        const std::complex<float> tmp01 = tmp0 + tmp1;
                        const std::complex<float> tmp23 = tmp2 + tmp3;
                        const std::complex<float> sum_field_terms = (tmp01)*(txft[0]) + (tmp23)*(txft[1]);
                        /*
                         Then, compute all the location terms combined. In 7.5-28, all the 
                         location terms are combined as below (excluding Doppler term):
                        
                                            T                            T
                                     /     r_rx*rx_dbar  \        /     r_tx*tx_dbar  \    
                                 exp| j2pi--------------- | * exp| j2pi--------------- |
                                    \       lambda_0    /        \       lambda_0    / 
                    
                        */
                        const std::complex<float> loc_terms = (*rxloc++)*(*txloc++);
                        *h++ = loc_terms*sum_field_terms*ksqrtinv[i]*norm_path_and_output; // repeating ksqrtinv is beneficial here

                        rxft += NUM_FIELD_SPHERICAL_ANGLES; // ptr will go out of bound by exactly 1 element, which is fine
                        txft += NUM_FIELD_SPHERICAL_ANGLES; // ptr will go out of bound by exactly 1 element, which is fine
                    }
                }
            }
        }
    } // namespace nlos

    namespace los {

        static void compute_cluster_gain(
            StaticPathGainConfig& los_gain_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const pdp::DelayProfileConfig& pdp_conf,
            const std::vector<float>& Phi
        ) {
            // Compute cluster gain for LoS component.

            const pdp::CDLCluster pdp_LoS_cluster = pdp_conf.Table[0];
            
            // 7.5 Eq. 7.5-13 for AoA and AoD (ray_offset = 0 for LoS)
            // Also wrap AoA and AoD to [-180, 180]
            const float phi_AoA = WRAP_AZIMUTH_ANGLE_180(pdp_LoS_cluster.aoa_deg);
            const float phi_AoD = WRAP_AZIMUTH_ANGLE_180(pdp_LoS_cluster.aod_deg);
            
            // 7.5 Eq. 7.5-18 for ZoA and ZoD (ray_offset = 0 for LoS)
            const float theta_ZoA = WRAP_ZENITH_ANGLE(pdp_LoS_cluster.zoa_deg);
            const float theta_ZoD = WRAP_ZENITH_ANGLE(pdp_LoS_cluster.zod_deg);
            
            // TODO: Angle scaling section 7.7.5.1 to be implemented

            // Get wavelength
            const float lambda_0 = ant_sys_conf.Wavelength();
            const float lambda_0_inv = 1/lambda_0;

            // -------------Processing for Tx----------------------------------------------------------------------------------------
            const antenna::geometry::los::SphericalDirection tx_sph_dir_conf = 
                antenna::geometry::los::get_spherical_unit_vector(phi_AoD, theta_ZoD); // angles of departure

            // Save Rx rhat for subsequent time-varying calculations
            los_gain_conf.TxSphericalUnitVector.assign(
                tx_sph_dir_conf.rhat.begin(),
                tx_sph_dir_conf.rhat.end()
            );

            // Get Tx antenna/subarray location vector
            const antenna::AntennaArrayConfig& tx_ant_conf = ant_sys_conf.TransmitAntennaArray;
            const std::vector<float>& tx_radiator_pos = tx_ant_conf.FieldPattern.ElementPositions;
            const antenna::geometry::PositionConfig& tx_pos_conf = tx_ant_conf.Position;

            // Calculate the location vector dbar of Tx
            std::vector<float> tx_dbar(tx_radiator_pos.size());
            for (std::size_t i = 0; i < tx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
                tx_dbar[i] = tx_radiator_pos[i]*lambda_0 + tx_pos_conf.x;
                tx_dbar[i + 1] = tx_radiator_pos[i + 1]*lambda_0 + tx_pos_conf.y;
                tx_dbar[i + 2] = tx_radiator_pos[i + 2]*lambda_0 + tx_pos_conf.z;
            }
            
            // Allocation for Tx field term and location term
            const std::size_t num_Tx = ant_sys_conf.NumInputSignals;
            std::vector<float> tx_field_term(num_Tx*NUM_FIELD_SPHERICAL_ANGLES);
            std::vector<std::complex<float>> tx_loc_term(num_Tx*NUM_RAYS_PER_CLUSTER);

            // Compute field and location terms
            for (std::size_t s = 0; s < num_Tx; ++s) {

                float* field_term = tx_field_term.data() + s*NUM_FIELD_SPHERICAL_ANGLES;

                antenna::field_pattern::los::compute_field_term(
                    field_term, 
                    tx_ant_conf, 
                    tx_sph_dir_conf, 
                    s
                );

                std::complex<float>* loc_term = tx_loc_term.data() + s*NUM_RAYS_PER_CLUSTER;

                antenna::field_pattern::los::get_location_term(
                    loc_term,
                    tx_sph_dir_conf.rhat,
                    tx_dbar,
                    lambda_0_inv,
                    s
                );
            }
            //------------End Processing for Tx----------------------------------------------------------------------------------------


            //--------------Processing for Rx------------------------------------------------------------------------------------------
            const antenna::geometry::los::SphericalDirection rx_sph_dir_conf = 
                antenna::geometry::los::get_spherical_unit_vector(phi_AoA, theta_ZoA); // angles of arrival

            // Save Rx rhat for subsequent time-varying calculations
            los_gain_conf.RxSphericalUnitVector.assign(
                rx_sph_dir_conf.rhat.begin(),
                rx_sph_dir_conf.rhat.end()
            );

            // Get Rx antenna/subarray location vector
            const antenna::AntennaArrayConfig&rx_ant_conf = ant_sys_conf.ReceiveAntennaArray;
            const std::vector<float>& rx_radiator_pos = rx_ant_conf.FieldPattern.ElementPositions;
            const antenna::geometry::PositionConfig& rx_pos_conf = rx_ant_conf.Position;

            // Calculate the location vector dbar of Rx
            std::vector<float> rx_dbar(rx_radiator_pos.size());
            for (std::size_t i = 0; i < rx_radiator_pos.size(); i += 3) { // +3 stride over pos(x,y,z)
                rx_dbar[i] = rx_radiator_pos[i]*lambda_0 + rx_pos_conf.x;
                rx_dbar[i + 1] = rx_radiator_pos[i + 1]*lambda_0 + rx_pos_conf.y;
                rx_dbar[i + 2] = rx_radiator_pos[i + 2]*lambda_0 + rx_pos_conf.z;
            }

            // Allocation for Tx field term and location term
            const std::size_t num_Rx = ant_sys_conf.NumOutputSignals;
            std::vector<float> rx_field_term(num_Rx*NUM_FIELD_SPHERICAL_ANGLES);
            std::vector<std::complex<float>> rx_loc_term(num_Rx*NUM_RAYS_PER_CLUSTER);

            // Compute field and location terms
            for (std::size_t u = 0; u < num_Rx; ++u) {

                float* field_term = rx_field_term.data() + u*NUM_FIELD_SPHERICAL_ANGLES;

                antenna::field_pattern::los::compute_field_term(
                    field_term, 
                    rx_ant_conf,
                    rx_sph_dir_conf, 
                    u
                );

                std::complex<float>* loc_term = rx_loc_term.data() + u*NUM_RAYS_PER_CLUSTER;
                
                antenna::field_pattern::los::get_location_term(
                    loc_term,
                    rx_sph_dir_conf.rhat,
                    rx_dbar,
                    lambda_0_inv,
                    u
                );
            }
            //------------End Processing for Rx----------------------------------------------------------------------------------------

            // Calculate polarization matrix
            const float XPR = pdp_conf.XPR;
            std::vector<std::complex<float>> polterm = calculate_polarization_matrix(
                Phi, 
                XPR, 
                pdp_conf.ClusterTypes, 
                PropagationCondition::LOS
            );

            // Calculate output normalization gain
            const float norm_output = channel_ctrl_conf.NormalizeChannelOutput ? 1.0f/std::sqrt(num_Rx) : 1.0f;

            // Get active ray and gain scaling
            // LoS active ray is always 1
            const float ksqrtinv = std::pow(10.0f, pdp_LoS_cluster.power_db*0.05f); // sqrt(10^pdB/10)

            // Calculate path normalization gain
            float norm_path_gain = 1.0f;

            if (channel_ctrl_conf.NormalizePathGains) {

                float avg_path_gains = 0.0f;

                for (const auto& p : pdp_conf.Table) {
                    avg_path_gains += std::pow(10.0f, p.power_db*0.1f);
                }

                norm_path_gain = 1.0f/std::sqrt(avg_path_gains);
            }
            
            // Calculate the MONSTROUS Eq. 7.5-28 excluding Doppler (the
            // last term). Doppler is to be calculated in time-varying channel response
            los_gain_conf.RayResponse.resize(num_Rx*num_Tx);
            std::complex<float>* h = los_gain_conf.RayResponse.data();
            // This implmentation doesn't include matrix multiplication as in 7.5-28 but an expanded
            // vectorized version for optimization.
            // For this LoS, we really don't need to iterate across the rays, since all rays have uniform equal power.
            // All it matters are the polarization matrix (polterm)

#if 0
            for (std::size_t s = 0; s < num_Tx; ++s) {

                const float* txft1 = tx_field_term.data() + s*NUM_FIELD_SPHERICAL_ANGLES; // F_theta
                const float* txft2 = txft1 + 1; // F_phi
                const std::complex<float>* txloc = tx_loc_term.data() + s*NUM_RAYS_PER_CLUSTER;

                for (std::size_t u = 0; u < num_Rx; ++u) {

                    const float* rxft1 = rx_field_term.data() + u*NUM_FIELD_SPHERICAL_ANGLES; // F_theta
                    const float* rxft2 = rxft1 + 1; // F_phi
                    const std::complex<float>* rxloc = rx_loc_term.data() + u*NUM_RAYS_PER_CLUSTER;
                    /*
                        (thth): theta-theta, (thph): theta-phi
                        ksqrtinv: 1/sqrt(kappa)

                    fieldTerms =

                            T
                    [F_rx_th]   [exp(jPhi_thth)          ksqrinv*exp(jPhi_thph]   [F_tx_th]
                    |       | * |                                             | * |       |
                    [F_rx_ph]   [ksqrinv*exp(jPhi_phth)         exp(jPhi_phph)]   [F_rx_ph]
                        
                    */

                    // This implmentation doesn't include matrix multiplication as in 7.5-28 
                    // but an expanded vectorized version for optimization.
                    // ksqrtinv included later
                    const std::complex<float> tmp0 = (*rxft1)*polterm[0];
                    const std::complex<float> tmp1 = (*rxft2)*polterm[1];
                    const std::complex<float> tmp2 = (*rxft1)*polterm[2];
                    const std::complex<float> tmp3 = (*rxft2)*polterm[3];
                    const std::complex<float> tmp01 = tmp0 + tmp1;
                    const std::complex<float> tmp23 = tmp2 + tmp3;
                    const std::complex<float> sum_field_terms = (tmp01)*(*txft1) + (tmp23)*(*txft2);

                    /*
                    Then, compute all the location terms combined. In 7.5-28, all the 
                    location terms are combined as below (excluding Doppler term)
                    
                                       T                            T
                                /     r_rx*rx_dbar  \        /     r_tx*tx_dbar  \    
                            exp| j2pi--------------- | * exp| j2pi--------------- |
                                \       lambda_0    /        \       lambda_0    / 
                
                    */
                    // No need to increment rxloc or txloc for LoS because it's the same for all ray
                    // but different for each antenna
                    const std::complex<float> loc_terms = (*rxloc)*(*txloc); 

                    *h++ = loc_terms*sum_field_terms*norm_path_gain*norm_output;
                }
            }
#else
            for (std::size_t u = 0; u < num_Rx; ++u) {

                const float* rxft1 = rx_field_term.data() + u*NUM_FIELD_SPHERICAL_ANGLES; // F_theta
                const float* rxft2 = rxft1 + 1; // F_phi
                const std::complex<float>* rxloc = rx_loc_term.data() + u*NUM_RAYS_PER_CLUSTER;

                const std::complex<float> tmp0 = (*rxft1)*polterm[0];
                const std::complex<float> tmp1 = (*rxft2)*polterm[1];
                const std::complex<float> tmp2 = (*rxft1)*polterm[2];
                const std::complex<float> tmp3 = (*rxft2)*polterm[3];
                const std::complex<float> tmp01 = tmp0 + tmp1;
                const std::complex<float> tmp23 = tmp2 + tmp3;
                const std::complex<float> rxloc_scaled = (*rxloc)*ksqrtinv;

                for (std::size_t s = 0; s < num_Tx; ++s) {

                    const float* txft1 = tx_field_term.data() + s*NUM_FIELD_SPHERICAL_ANGLES; // F_theta
                    const float* txft2 = txft1 + 1; // F_phi
                    const std::complex<float>* txloc = tx_loc_term.data() + s*NUM_RAYS_PER_CLUSTER;
                    // First, compute all the field terms combined. In 7.5-28, all the 
                    // field terms are combined as below:
                    /*
                        (thth): theta-theta, (thph): theta-phi
                        ksqrtinv: 1/sqrt(kappa)

                     fieldTerms =

                              T
                     [F_rx_th]   [exp(jPhi_thth)          ksqrinv*exp(jPhi_thph]   [F_tx_th]
                     |       | * |                                             | * |       |
                     [F_rx_ph]   [ksqrinv*exp(jPhi_phth)         exp(jPhi_phph)]   [F_rx_ph]
                        
                    */

                    // This implmentation doesn't include matrix multiplication as in 7.5-28 
                    // but an expanded vectorized version for optimization.
                    // ksqrtinv included in rxloc_scaled
                    const std::complex<float> sum_field_terms = (tmp01)*(*txft1) + (tmp23)*(*txft2);

                    /*
                     Then, compute all the location terms combined. In 7.5-28, all the 
                     location terms are combined as below (excluding Doppler term):
                    
                                        T                            T
                                 /     r_rx*rx_dbar  \        /     r_tx*tx_dbar  \    
                             exp| j2pi--------------- | * exp| j2pi--------------- |
                                 \       lambda_0    /        \       lambda_0    / 
                
                    */
                    // No need to increment rxloc or txloc for LoS because it's the same for all ray
                    // but different for each antenna
                    *h++ = rxloc_scaled*(*txloc)*sum_field_terms*norm_path_gain*norm_output;
                }
            }
#endif
        }
    } // namespace los

    static std::vector<std::complex<float>> calculate_polarization_matrix(
        const std::vector<float>& Phi,
        const float XPR,
        const std::vector<PropagationCondition> all_cluster_type,
        const PropagationCondition this_cluster_type
    ) {
        // Calculate polarization matrix for the specified cluster type
        /* 
         The Phi (Angular interaction) is arranged like this:
         | theta-theta |
         |   phi-theta |
         | theta-phi   |
         |   phi-phi   |
        */
        
        std::size_t num_cluster = all_cluster_type.size();
        std::size_t this_cluster_count = 0;
        for (std::size_t i = 0; i < num_cluster; ++i) {
            if (all_cluster_type[i] == this_cluster_type) {
                ++this_cluster_count;
            }
        }

        std::vector<std::complex<float>> polmat;
        std::complex<float>* polmat_ptr;

        switch (this_cluster_type) {

            case PropagationCondition::LOS: {

                // Reallocate polmat and update its pointer
                polmat.resize(NUM_PATH_ANGLES);
                polmat_ptr = polmat.data();

                // LoS is always at the first index if available
                std::complex<float> cphi = std::complex(
                    std::cos(Phi[0]),
                    std::sin(Phi[0])
                );

                // For LoS case, it is really simple
                // only theta-theta and phi-phi matter, no cross angle
                polmat_ptr[0] =  cphi;
                polmat_ptr[3] = -cphi; // phi-phi is negative
                break;
            }

            case PropagationCondition::NLOS: {

                // Number of path angle samples along all NLoS clusters along all rays
                const std::size_t num_scr = this_cluster_count*NUM_RAYS_PER_CLUSTER*NUM_PATH_ANGLES;

                // Reallocate polmat and update its pointer
                polmat.resize(num_scr);
                polmat_ptr = polmat.data();

                // Number of samples along all clusters including any LoS along all rays
                const std::size_t num_samp = num_cluster*NUM_RAYS_PER_CLUSTER;

                // Compute kappa: the cross-polarization power ratio (XPR), given in dB.
                // XPR is defined as the ratio of power in co-polarization (theta-theta or phi-phi)
                // to cross-polarization (theta-phi or phi-theta). See TR 38.901 Eq 7.5-21.
                //sqrt of 1/kappa combined
                const float k = std::pow(10.0f, -XPR*0.05f); // 1/sqrt(10^(XPR/10))

                // If LoS exist, begin_idx will start with 1 index later 
                // because index 0 will always be LoS
                const bool hasLoS = this_cluster_count == num_cluster - 1;
                const std::size_t begin_idx = static_cast<std::size_t>(hasLoS);

                // polarization pair stride (index bookkeeping) for initial phase Phi flat array
                const std::size_t s1 = num_cluster*NUM_RAYS_PER_CLUSTER;;
                const std::size_t s2 = s1*2;
                const std::size_t s3 = s1*3;

                // Calculate complex exponentials of initial phases.
                // Phi is expected to be 4*N (for theta-theta, theta-phi, phi-theta, phi-phi).
                // Outer loop iterating over all clusters including any LoS, 
                // while incrementing conditionally based on LoS existence
                for (std::size_t m = 0; m < num_samp; m += (this_cluster_count + begin_idx)) {
                    
                    // Inner loops begins with the outer loop value as offset 
                    // and iterating across number of NLoS
                    for (std::size_t i = m; i < this_cluster_count + m; ++i) {

                        const std::size_t im = i + begin_idx;
                        
                        // Calculate complex exponentials of initial phases.
                        // Phi is expected to be 4×N (for theta-theta, theta-phi, phi-theta, phi-phi).
                        // For polarization matrix, we apply per-ray scaling based on kappa.
                        // See Eq 7.5-22 Phi
                        // Cross-pol tp and pt will be scaled by sqrt_kappa_in
                        const std::complex<float> tt = {std::cos(Phi[im]), std::sin(Phi[im])};
                        const std::complex<float> tp = {k*std::cos(Phi[im + s1]), k*std::sin(Phi[im + s1])};
                        const std::complex<float> pt = {k*std::cos(Phi[im + s2]), k*std::sin(Phi[im + s2])};
                        const std::complex<float> pp = {std::cos(Phi[im + s3]), std::sin(Phi[im + s3])};

                        *polmat_ptr++ = tt; // theta-theta: Phi^(theta,theta)_(n,m)
                        *polmat_ptr++ = pt; // phi-theta:   Phi^(phi,theta)_(n,m)
                        *polmat_ptr++ = tp; // theta-phi:   Phi^(theta,phi)_(n,m)
                        *polmat_ptr++ = pp; // phi-phi:     Phi^(phi,phi)_(n,m)
                    }
                }
                break;
                
            }
            case channels::PropagationCondition::SUBCLUSTER_NLOS:

                
            default:

        }
        return polmat;
    }

} //namespace channels::cdl::response
