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

    /***************** Function Declarations ********************/
    static std::vector<float> get_sample_times(
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    static std::vector<std::complex<float>> compute_doppler_term(
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
        // Generating time-varying path gains by adding doppler effects to static path gains

        tv_gain_conf.TimeVector = get_sample_times(channel_ctrl_conf, ant_sys_conf, mobility_conf);

        const std::uint16_t num_input = ant_sys_conf.NumInputSignals;
        const bool HasLoS = pdp_conf.HasLoS;

        if (HasLoS) {
            const std::vector<float>& tx_los_rhat = los_gain_conf.TxSphericalUnitVector;
            const std::vector<float>& rx_los_rhat = los_gain_conf.RxSphericalUnitVector;
            const std::vector<std::complex<float>>& los_H_static = los_gain_conf.Response;
        }

        const std::vector<float>& tx_nlos_rhat = nlos_gain_conf.TxSphericalUnitVector;
        const std::vector<float>& rx_nlos_rhat = nlos_gain_conf.RxSphericalUnitVector;
        const std::vector<std::complex<float>>& nlos_H_static = nlos_gain_conf.Response;








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

    static std::vector<std::complex<float>> compute_doppler_term(
        const bool HasLoS,
        const std::vector<float>& time_vect,
        const StaticPathGainConfig& los_gain_conf,
        const StaticPathGainConfig& nlos_gain_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    ) {
        // Compute Doppler term in equation 7.5-28/29 in step 11
        /*
                        T
           /        rhat_(rx)_(LOS/NLOS) * vbar    \
        exp|j2*pi--------------------------------*t|
           \                 lambda_0              /

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
        if (ut_dot.size() == 2) {

            rx_theta_v = tx_theta_v = ut_dot[1]; // azimuth
            rx_phi_v = tx_phi_v = 0.0f; // zenith

        } else if (ut_dot.size() == 4) {

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

        std::vector<float> tx_vbar, rx_vbar;
        for (std::size_t i = 0; i < antenna::geometry::NUM_3D_AXIS; ++i) {
            tx_vbar[i] = tx_v*tx_v_sph_dir.rhat[i];
            rx_vbar[i] = rx_v*rx_v_sph_dir.rhat[i];
        }

        const mobility::ScattererConfig& scatterer_conf = mobility_conf.Scatterer;
        const std::size_t num_scatterer = scatterer_conf.States.size();
        const std::vector<float>& speed = scatterer_conf.Speed;
        const std::vector<bool>& states = scatterer_conf.States;

        float tx_los_rhat_x = 0.0f;
        float tx_los_rhat_y = 0.0f;
        float tx_los_rhat_z = 0.0f;

        float rx_los_rhat_x = 0.0f;
        float rx_los_rhat_y = 0.0f;
        float rx_los_rhat_z = 0.0f;

        if (HasLoS) {
            float tx_los_rhat_x = los_gain_conf.TxSphericalUnitVector[0];
            float tx_los_rhat_y = los_gain_conf.TxSphericalUnitVector[1];
            float tx_los_rhat_z = los_gain_conf.TxSphericalUnitVector[2];

            float rx_los_rhat_x = los_gain_conf.RxSphericalUnitVector[0];
            float rx_los_rhat_y = los_gain_conf.RxSphericalUnitVector[1];
            float rx_los_rhat_z = los_gain_conf.RxSphericalUnitVector[2];
        
        } 

        std::size_t num_time_pts = time_vect.size();
        std::vector<std::complex<float>> tv_grid(num_time_pts*num_scatterer);
        std::complex<float>* tv = tv_grid.data();

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

            // Calculate tx_vnm
            const float tx_vnm_x = (*tx_rh_nlos_ptr++)*tx_vbar[0] + tx_los_rhat_x;
            const float tx_vnm_y = (*tx_rh_nlos_ptr++)*tx_vbar[1] + tx_los_rhat_y;
            const float tx_vnm_z = (*tx_rh_nlos_ptr++)*tx_vbar[2] + tx_los_rhat_z;
            const float tx_vnm   = tx_vnm_x + tx_vnm_y + tx_vnm_z;

            // Calculate rx_vnm
            const float rx_vnm_x = (*rx_rh_nlos_ptr++)*rx_vbar[0] + rx_los_rhat_x;
            const float rx_vnm_y = (*rx_rh_nlos_ptr++)*rx_vbar[1] + rx_los_rhat_y;
            const float rx_vnm_z = (*rx_rh_nlos_ptr++)*rx_vbar[2] + rx_los_rhat_z;
            const float rx_vnm   = rx_vnm_x + rx_vnm_y + rx_vnm_z;

            const float tmp = (tx_rx_scatterer + tx_vnm + rx_vnm)*lambda_0_inv;
            const float tmp1 = TWO_PI*tmp;

            for (const auto& t : time_vect) {

                const float tmp2 = tmp1*t;

                *tv++ = {std::cos(tmp2), std::sin(tmp2)};
            }
        }
        return tv_grid;
    }

} // namespace channels::cdl::response