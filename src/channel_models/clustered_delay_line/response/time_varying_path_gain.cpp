#include <variant>
#include <vector>
#include <cmath>
#include <variant>
#include "generator.hpp"
#include "../mobility/mobility.hpp"


namespace channels::cdl::response {

    static std::vector<float> get_sample_times(
        const std::uint64_t num_time_samples,
        const std::uint16_t num_input,
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    void generate_time_varying_path_gains(
            TimeVaryingPathGainConfig& tv_gain_conf,
            const StaticPathGainConfig& los_gain_conf,
            const StaticPathGainConfig& nlos_gain_conf,
            const ChannelControlConfig& channel_ctrl_conf,
            const antenna::AntennaSystemConfig& ant_sys_conf,
            const mobility::MobilityConfig& mobility_conf
    ) {
        // Generating time-varying path gains by adding doppler effects to static path gains

        const std::uint64_t samp_size = channel_ctrl_conf.NumTimeSamples;
        const std::uint16_t num_input = ant_sys_conf.NumInputSignals;





    }

    static std::vector<float> get_sample_times(
        const std::uint64_t num_time_samples,
        const std::uint16_t num_input,
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    ) {
        // Return time samples based on channel control config and
        // cofficient generation sampling rate F_cg

        const float fc = ant_sys_conf.CarrierFrequency;
        const float time_offset = channel_ctrl_conf.CurrentTime;
        std::variant<float, std::vector<float>> max_doppler;

        const bool is_max_doppler_scalar = std::holds_alternative<float>(mobility_conf.MaximumDopplerShift);
        
        if (is_max_doppler_scalar) {
            max_doppler = std::get<float>(mobility_conf.MaximumDopplerShift);
        } else {
            max_doppler = std::get<std::vector<float>>(mobility_conf.MaximumDopplerShift);
        }





    }

}