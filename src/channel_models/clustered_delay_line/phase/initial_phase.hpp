#pragma once

#include "../antenna_array/antenna.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"

namespace channels::cdl::phase {

    std::vector<float> generateInitialPhase(
        pdp::DelayProfileConfig& pdp_conf,
        RandomStreamConfig& rstream_conf,
        const ChannelControlConfig& control_conf,
        const antenna::AntennaSystemConfig& ant_conf
    );

}