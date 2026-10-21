#include <cstdint>
#include <ranges>
#include <variant>
#include <vector>
#include <array>
#include <cmath>
#include <variant>
#include "generator.hpp"
#include "../mobility/mobility.hpp"

/************************* Constants *****************************/
static constexpr float DOUBLE_BOUNCE_PATH = 2.0f;
static constexpr float TWO_PI = M_PI*2.0f;

namespace channels::cdl::response {

    /***************** Local Type Definitions ********************/
    struct DopplerTerm {
        std::vector<std::complex<float>> los_term;
        std::vector<std::complex<float>> nlos_term;
    };

    /***************** Function Declarations ********************/
    static std::vector<float> get_sample_times(
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    static DopplerTerm compute_doppler_term(
        const bool HasLoS,
        const std::vector<float>& time_vect,
        const StaticPathGainConfig& los_gain_conf,
        const StaticPathGainConfig& nlos_gain_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    /***************** Function Implementations ********************/
    void generate_time_varying_path_gains(
            TimeVaryingPathGainConfig& tv_gain_conf,
            const StaticPathGainConfig& los_gain_conf,
            const StaticPathGainConfig& nlos_gain_conf,
            const pdp::DelayProfileConfig& pdp_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const mobility::MobilityConfig& mobility_conf
    ) {
        // Generating time-varying path gains by adding doppler effects to static path gains.

        tv_gain_conf.TimeVector = get_sample_times(channel_ctrl_conf, ant_sys_conf, mobility_conf);
        const std::size_t num_samp_times = tv_gain_conf.TimeVector.size();

        const uint8_t has_los = static_cast<uint8_t>(pdp_conf.HasLoS);
        const std::vector<std::complex<float>>& los_H_static = los_gain_conf.Response;
        const std::vector<std::complex<float>>& nlos_H_static = nlos_gain_conf.Response;

        const DopplerTerm doppler_conf = compute_doppler_term(
            has_los,
            tv_gain_conf.TimeVector,
            los_gain_conf,
            nlos_gain_conf,
            ant_sys_conf,
            mobility_conf
        );

        const std::size_t numRx = ant_sys_conf.NumOutputSignals;
        const std::size_t numTx = ant_sys_conf.NumInputSignals;
        const std::size_t num_all_clusters = pdp_conf.ClusterTypes.size();
        const std::size_t num_all_samples = num_all_clusters*NUM_RAYS_PER_CLUSTER;
        const std::size_t num_nlos_clusters = num_all_clusters - has_los;
        const std::size_t num_nlos_samples = num_all_samples - NUM_RAYS_PER_CLUSTER;
        const std::size_t num_nlos_gain_samples = num_nlos_samples*numRx*numTx;
        const std::complex<float>* nlos_dp = doppler_conf.nlos_term.data();
        const std::vector<std::complex<float>>& los_dp_term = doppler_conf.los_term;
        const std::vector<std::complex<float>>& nlos_dp_term = doppler_conf.nlos_term;
        

        std::vector<std::complex<float>> H(num_samp_times*num_all_clusters*numTx*numRx);
        std::complex<float>* h = H.data();
        // Minus num_samp_times to reset the h ptr back to the 0th sample for next stride
        const std::size_t los_h_stride = num_all_clusters*num_samp_times - num_samp_times;

        std::vector<std::complex<float>> H2(num_samp_times*num_all_clusters*NUM_RAYS_PER_CLUSTER*numTx*numRx);
        std::complex<float>* h2 = H2.data();

        // Combine static and doppler term to produce time-varying channel
        if (pdp_conf.HasLoS) {

            // LoS component processing
            for (std::size_t u = 0; u < numRx; ++u ) {

                for (std::size_t s = 0; s < numTx; ++s) {
                    const std::complex<float> gain = los_H_static[u*numTx + s];

                    for (std::size_t t = 0; t < num_samp_times; ++t) {
                        *h++ = gain*los_dp_term[t];
                    }
                    h += los_h_stride;
                }
            }

            // reset h ptr to beginning of nlos samples
            std::vector<std::size_t> ss;
            h = H.data() + num_samp_times;
            std::vector<std::complex<float>> ray_accum(num_samp_times);

            for (std::size_t u = 0; u < numRx; ++u ) {

                std::size_t u_stride = u*num_nlos_samples*numTx;
                for (std::size_t s = 0; s < numTx; ++s) {

                    std::size_t s_stride = s*num_nlos_samples;
                    for (std::size_t i = 0; i < num_nlos_clusters; ++i) {

                        for (std::size_t m = 0; m < NUM_RAYS_PER_CLUSTER; ++m) {
                            std::size_t m_stride = m*num_nlos_clusters;
                            const std::size_t doppler_base_idx = num_nlos_clusters*num_samp_times*m + i*num_samp_times;
                            std::size_t sg_base_idx_2 = u_stride + s_stride + m_stride + i;
                            const std::complex<float> gain = nlos_H_static[sg_base_idx_2];

                            for (std::size_t t = 0; t < num_samp_times; ++t) {
                                const std::complex<float> C = gain*nlos_dp_term[doppler_base_idx + t];
                                ray_accum[t] += C;
                            }
                        }

                        // Load the accumulated rays into H array
                        for (std::size_t t = 0; t < num_samp_times; ++t) {
                            *h++ = ray_accum[t];
                            ray_accum[t] = {0.0f, 0.0f}; // reset ray accumlation
                        }
                    }
                    h += num_samp_times; // offset ptr by num_samp_times

/*
                    for (std::size_t m = 0; m < NUM_RAYS_PER_CLUSTER; ++m) {

                        std::complex<float> ray_accum = 0;

                        for (std::size_t i = 0; i < num_nlos_clusters; ++i) {
                            const std::size_t cluster_per_ray = m*num_nlos_clusters;
                            const std::size_t doppler_base_idx = i*num_samp_times + cluster_per_ray;
                            const std::size_t sg_base_idx = u*numTx + s + cluster_per_ray + i;
                            const std::complex<float> gain = nlos_H_static[sg_base_idx];

                            for (std::size_t t = 0; t < num_samp_times; ++t) {

                                *h++ += gain*nlos_dp_term[doppler_base_idx + t];

                            }
                        }

                    }
                        */
                }
            }

            std::cout << "Done!!!" << std::endl;

        } else {

        }

    }

    static std::vector<float> get_sample_times(
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    ) {
        // Return time samples based on channel control config and
        // cofficient generation sampling rate F_cg.

        const mobility::ScattererConfig& scatterer_conf = mobility_conf.Scatterer;
        const float fc = ant_sys_conf.CarrierFrequency;
        const float time_offset = channel_ctrl_conf.CurrentTime;
        const bool is_max_doppler_scalar = std::holds_alternative<float>(mobility_conf.MaximumDopplerShift);

        float tx_rx_doppler;
        float max_scatter_v;
        if (is_max_doppler_scalar) {
            // If input is a scalar float, then it's the maximum Doppler shift of the Rx.
            tx_rx_doppler = std::get<float>(mobility_conf.MaximumDopplerShift);
            max_scatter_v = 0.0f;
        } else {
            // If inpyt is a 2-element array then sum the maximum Doppler shift of the Rx and T.
            const auto& v = std::get<std::array<float, 2>>(mobility_conf.MaximumDopplerShift);
            tx_rx_doppler = v[0] + v[1];
            max_scatter_v = scatterer_conf.MaximumScattererSpeed;
        }

        const float fs = channel_ctrl_conf.SampleRate;
        const float samp_density = channel_ctrl_conf.SampleDensity;
        const std::uint64_t num_wf_samp = channel_ctrl_conf.NumTimeSamples;
        const float lambda_0 = ant_sys_conf.Wavelength();

        // Doppler from moving scatterer, causing both the forward and reflected
        // path, like a double bounce.
        const float scatterer_doppler = DOUBLE_BOUNCE_PATH*max_scatter_v/lambda_0; // Hz

        // The number 2 here indicates F_cg must sample at least 2 times the max
        // Doppler freq as per Nyquist theorem but can be oversampled by the sampDensity factor.
        const float osf = 2*channel_ctrl_conf.SampleDensity;

        // Then the coefficient generation sampling rate F_cg is:
        float F_cg = osf*(tx_rx_doppler + scatterer_doppler); // samps/s

        // Calculate total waveform duration in seconds:
        const float time_wf_sec = num_wf_samp/fs;

        // Calculate the number of channel coefficient samples needed to cover
        // the waveform duration, sampled at F_cg (channel generation rate).
        // Pad by at least half a sample period to ensure ZOH has one final sample to hold
        // Thisprevents end-of-waveform interpolation artifacts (important for channel filter later).
        uint64_t num_time_pts = static_cast<uint64_t>(std::ceil(time_wf_sec*F_cg + 0.5f));

        // Oversample channel taps to resolve Doppler effects (2*fD*SampleDensity),
        // but clip to SampleRate since that's the max resolution the waveform
        // allows. Physics first, practical constraints second.
        F_cg = std::min(F_cg, fs);
        const float F_cg_inv = 1.0f/F_cg;
        num_time_pts = std::min(num_time_pts, num_wf_samp);

        // Calculate time vectors at F_cg
        std::vector<float> time_vect(num_time_pts);
        if (num_time_pts > 1) {
            for (std::size_t i = 0; i < num_time_pts; ++i) {
                time_vect[i] = i*F_cg_inv + time_offset;
            }
        }
        return time_vect;
    }

    static DopplerTerm compute_doppler_term(
        const bool has_los,
        const std::vector<float>& time_vect,
        const StaticPathGainConfig& los_gain_conf,
        const StaticPathGainConfig& nlos_gain_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    ) {
        // Compute the Doppler term in Step 11, Eq. (7.5-28)/(7.5-29) of 3GPP TR 38.901.
        /*
         The original equations assume a stationary transmitter (Tx) and
         account only for receiver (Rx) mobility:
        
                             T
                 /       rhat_(rx)_(LOS/NLOS) * vbar_(rx)  \
             exp|j2*pi*----------------------------------*t|
                 \                 lambda_0                /
        
         However, when both Tx and Rx are mobil (dual mobility), the Doppler shift must
         include the contributions from both ends of the radio link,
         as specified in Section 7.6.10 (Dual Mobility). 
         This implementation supports dual mobility scenarios.
        
         The resulting Doppler term is:
        
                              T                                  T
                 /        rhat_(rx)_(LOS/NLOS) * vbar_(rx) + rhat_(tx)_(LOS/NLOS) * vbar_(tx)    \
             exp|j2*pi*------------------------------------------------------------------------*t|
                 \                                    lambda_0                                   /
        
         where:
           - rhat_(rx): Unit vector corresponding to the ray's arrival direction (AoA/ZoA) at the receiver.
           - rhat_(tx): Unit vector corresponding to the ray's departure direction (AoD/ZoD) at the transmitter.
           - vbar_(rx): Receiver velocity vector.
           - vbar_(tx): Transmitter velocity vector.
           - lambda_0: Carrier wavelength.
           - t:        Time instant.
        
         Each velocity is projected onto its corresponding ray direction
         to obtain the Doppler contribution from that side of the link.
         The two contributions are added in the exponent, following the
         3GPP departure/arrival direction conventions.
        
         For a stationary Tx (vbar_(tx) = 0), this expression reduces
         to the original Doppler term in Eq. (7.5-28)/(7.5-29).
        */

        const float lambda_0 = ant_sys_conf.Wavelength();
        const float lambda_0_inv = 1.0f/lambda_0;

        // Rx and Tx speeds
        const bool is_max_doppler_scalar = std::holds_alternative<float>(mobility_conf.MaximumDopplerShift);
        float tx_v, rx_v;

        if (is_max_doppler_scalar) {
            rx_v = std::get<float>(mobility_conf.MaximumDopplerShift)*lambda_0;
            tx_v = 0.0f;
        } else {
            const auto& v = std::get<std::array<float, 2>>(mobility_conf.MaximumDopplerShift);
            rx_v = v[0]*lambda_0;
            tx_v = v[1]*lambda_0;
        }

        // User-Terminal Direction of travel
        const std::vector<float>& ut_dot = mobility_conf.UTDirectionOfTravel;
        float rx_theta_v, rx_phi_v, tx_theta_v, tx_phi_v;
        rx_theta_v = tx_theta_v = ut_dot[1]; // azimuth
        rx_phi_v = tx_phi_v = 0.0f; // zenith

        if (ut_dot.size() == 4) {
            rx_theta_v = ut_dot[0];
            tx_theta_v = ut_dot[1];
            rx_phi_v = ut_dot[2];
            tx_phi_v = ut_dot[3];
        }

        // Use the los namespace get rhat function because it accepts scalar phi and theta
        const antenna::geometry::los::SphericalDirection tx_v_sph_dir = 
            antenna::geometry::los::get_spherical_unit_vector(tx_phi_v, tx_theta_v);

        const antenna::geometry::los::SphericalDirection rx_v_sph_dir = 
            antenna::geometry::los::get_spherical_unit_vector(rx_phi_v, rx_theta_v);

        std::vector<float> tx_vbar(antenna::geometry::NUM_3D_AXIS);
        std::vector<float> rx_vbar(antenna::geometry::NUM_3D_AXIS);
        for (std::size_t i = 0; i < antenna::geometry::NUM_3D_AXIS; ++i) {
            tx_vbar[i] = tx_v*tx_v_sph_dir.rhat[i];
            rx_vbar[i] = rx_v*rx_v_sph_dir.rhat[i];
        }

        // Only get the nlos scatterers because los cluster has no scatterer 
        // and is constant across all rays
        const mobility::ScattererConfig& scatterer_conf = mobility_conf.Scatterer;
        const std::size_t num_scatterer = scatterer_conf.NLoSStates.size();
        const std::vector<float>& speed = scatterer_conf.NLoSSpeed;
        const std::vector<uint8_t>& states = scatterer_conf.NLoSStates;

        DopplerTerm doppler_term_conf;
        std::size_t num_time_pts = time_vect.size();
        std::vector<std::complex<float>>& los_term = doppler_term_conf.los_term;
        los_term.resize(num_time_pts); // only one scatterer for LoS
        doppler_term_conf.nlos_term.resize(num_time_pts*num_scatterer);
        std::complex<float>* nlos_dp = doppler_term_conf.nlos_term.data();

        // Handle LoS doppler term
        float tx_los_rhat_x = 0.0f;
        float tx_los_rhat_y = 0.0f;
        float tx_los_rhat_z = 0.0f;
        float tx_los_rhat_vbar_sum = 0.0f;

        float rx_los_rhat_x = 0.0f;
        float rx_los_rhat_y = 0.0f;
        float rx_los_rhat_z = 0.0f;
        float rx_los_rhat_vbar_sum = 0.0f;

        if (has_los) {
            tx_los_rhat_x = los_gain_conf.TxSphericalUnitVector[0];
            tx_los_rhat_y = los_gain_conf.TxSphericalUnitVector[1];
            tx_los_rhat_z = los_gain_conf.TxSphericalUnitVector[2];
            tx_los_rhat_vbar_sum = tx_los_rhat_x*tx_vbar[0] + 
                                tx_los_rhat_y*tx_vbar[1] + 
                                tx_los_rhat_z*tx_vbar[2];

            rx_los_rhat_x = los_gain_conf.RxSphericalUnitVector[0];
            rx_los_rhat_y = los_gain_conf.RxSphericalUnitVector[1];
            rx_los_rhat_z = los_gain_conf.RxSphericalUnitVector[2];
            rx_los_rhat_vbar_sum = rx_los_rhat_x*rx_vbar[0] + 
                                rx_los_rhat_y*rx_vbar[1] + 
                                rx_los_rhat_z*rx_vbar[2];
            
            const float tmp0 = TWO_PI*(tx_los_rhat_vbar_sum + rx_los_rhat_vbar_sum)*lambda_0_inv; // no scatterer

            for (std::size_t i = 0; i < num_time_pts; ++i) {
                const float tmp1 = tmp0*time_vect[i];
                los_term[i] = {std::cos(tmp1), std::sin(tmp1)};
            }
        } 

        // Compute the doppler or doppler phase term
        const std::vector<float>& tx_nlos_rhat = nlos_gain_conf.TxSphericalUnitVector;
        const float* tx_rh_nlos_ptr = tx_nlos_rhat.data();
        const std::vector<float>& rx_nlos_rhat = nlos_gain_conf.RxSphericalUnitVector;
        const float* rx_rh_nlos_ptr = rx_nlos_rhat.data();
        const float lamb_inv_2pi = lambda_0_inv*TWO_PI;

        for (std::size_t i = 0; i < num_scatterer; ++i) {
            const auto& sp = speed[i]; // alias to avoid longer lines
            const auto& st = states[i];

            // Double bounce path (2) because Tx->Scatterer->Rx bounce
            const float tx_rx_scatterer = DOUBLE_BOUNCE_PATH*sp*static_cast<float>(st);

            // Calculate tx_vnm Eq. 7.6-45 (LoS) & 7.6-46 (NLoS)
            const float tx_vnm_x = (*tx_rh_nlos_ptr++)*tx_vbar[0];
            const float tx_vnm_y = (*tx_rh_nlos_ptr++)*tx_vbar[1];
            const float tx_vnm_z = (*tx_rh_nlos_ptr++)*tx_vbar[2];
            const float tx_vnm   = tx_vnm_x + tx_vnm_y + tx_vnm_z;

            // Calculate rx_vnm
            const float rx_vnm_x = (*rx_rh_nlos_ptr++)*rx_vbar[0];
            const float rx_vnm_y = (*rx_rh_nlos_ptr++)*rx_vbar[1];
            const float rx_vnm_z = (*rx_rh_nlos_ptr++)*rx_vbar[2];
            const float rx_vnm   = rx_vnm_x + rx_vnm_y + rx_vnm_z;

            const float tmp = (tx_rx_scatterer + tx_vnm + rx_vnm)*lambda_0_inv;
            const float tmp1 = TWO_PI*tmp;

            for (const auto& t : time_vect) {
                const float tmp2 = tmp1*t;
                *nlos_dp++ = {std::cos(tmp2), std::sin(tmp2)};
            }
        }
        return doppler_term_conf;
    }

} // namespace channels::cdl::response