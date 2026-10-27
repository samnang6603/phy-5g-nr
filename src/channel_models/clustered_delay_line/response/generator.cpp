#include <complex>
#include <cstddef>
#include "generator.hpp"

namespace channels::cdl::response {
 
    /***************** Function Implementations ********************/

    static void generate_path_gains(
        PathGainConfig& path_gain_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& Phi
    );

    //static void generate_ofdm_response();

    std::vector<std::complex<float>>& generate_response_output(
        PathGainConfig& path_gain_conf,
        OFDMResponseConfig& ofdm_response_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    ) {

        // Generate time-varying path gains
        generate_time_varying_path_gains(
            path_gain_conf.TimeVaryingGain,
            path_gain_conf.LoSStaticGain,
            path_gain_conf.NLoSStaticGain,
            pdp_conf,
            channel_ctrl_conf,
            ant_sys_conf,
            mobility_conf
        );

        switch (channel_ctrl_conf.ChannelResponseOutput) {

            case ChannelResponseOutputType::PATH_GAINS: {

                return path_gain_conf.TimeVaryingGain.ClusterResponse;
            }

            case ChannelResponseOutputType::OFDM_RESPONSE: {

                //throw std::logic_error(
                //    "OFDM response generation is not implemented."
                //); 
                generate_ofdm_response(
                    ofdm_response_conf,
                    path_gain_conf.TimeVaryingGain,
                    pdp_conf,
                    ant_sys_conf,
                    channel_ctrl_conf
                );
                
            }

            default: {

                throw std::invalid_argument(
                    "Unsupported channel response output type."
                );
            }
        }
    }

    static void generate_path_gains(
        PathGainConfig& path_gain_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& Phi
    ) {

        generate_static_path_gains(
            path_gain_conf.LoSStaticGain,
            path_gain_conf.NLoSStaticGain, 
            channel_ctrl_conf, 
            pdp_conf, 
            ant_sys_conf, 
            coupling, 
            Phi
        );

        generate_time_varying_path_gains(
            path_gain_conf.TimeVaryingGain, 
            path_gain_conf.LoSStaticGain, 
            path_gain_conf.NLoSStaticGain, 
            pdp_conf, 
            channel_ctrl_conf, 
            ant_sys_conf, 
            mobility_conf
        );
    }

} // namespace channels::cdl::response