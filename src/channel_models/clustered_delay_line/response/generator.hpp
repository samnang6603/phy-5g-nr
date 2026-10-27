#pragma once

#include <array>
#include "../antenna_array/antenna.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"
#include "../mobility/mobility.hpp"
#include "../../../modules/carrier/carrier.hpp"


// Ray offset in 7.5 Table 7.5-3
constexpr std::array<float, NUM_RAYS_PER_CLUSTER> RAY_OFFSET_ALPHA = {
    0.0447f,  -0.0447f,   0.1413f,  -0.1413f,   0.2492f,  
   -0.2492f,   0.3715f,  -0.3715f,   0.5129f,  -0.5129f,
    0.6797f,  -0.6797f,   0.8844f,  -0.8844f,  1.1481f,
   -1.1481f,   1.5195f,  -1.5195f,   2.1551f,  -2.1551f
};

namespace channels::cdl::response {

    struct StaticPathGainConfig {
        std::vector<std::complex<float>> RayResponse;
        std::vector<float> TxSphericalUnitVector;
        std::vector<float> RxSphericalUnitVector;
    };

    struct TimeVaryingPathGainConfig {
        std::vector<std::complex<float>> ClusterResponse;
        std::vector<float> TimeVector;
    };

    struct PathGainConfig {
        StaticPathGainConfig LoSStaticGain;
        StaticPathGainConfig NLoSStaticGain;
        TimeVaryingPathGainConfig TimeVaryingGain;
        std::vector<float> TimeSamples; // channel snapshot
    };

    struct OFDMResponseConfig {
        std::vector<std::complex<float>> Response;
        modules::Carrier Carrier;
    };

    void generate_static_path_gains(
        StaticPathGainConfig& los_gain_conf,
        StaticPathGainConfig& nlos_gain_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& Phi
    );
    
    void generate_time_varying_path_gains(
        TimeVaryingPathGainConfig& tv_gain_conf,
        const StaticPathGainConfig& los_gain_conf,
        const StaticPathGainConfig& nlos_gain_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    void generate_ofdm_response(
        OFDMResponseConfig& ofdm_response_conf,
        const TimeVaryingPathGainConfig& time_varying_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const ChannelControlConfig& channel_ctrl_conf
    );

    std::vector<std::complex<float>>& generate_response_output(
        PathGainConfig& path_gain_conf,
        OFDMResponseConfig& ofdm_response_conf,
        const ChannelControlConfig& channel_ctrl_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const mobility::MobilityConfig& mobility_conf
    );

    


    
}