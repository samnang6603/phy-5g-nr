#pragma once

#include "../cdl.hpp"
#include "../antenna_structure/antenna_structure.hpp"

namespace channels::cdl::phase {

    std::vector<float> generateInitialPhase(
        pdp::DelayProfileConfig& pdp_conf,
        nrCDLChannel::RandomStreamConfig& rstream_conf,
        nrCDLChannel::ChannelControlConfig& control_conf,
        antenna::AntennaSystemConfig& ant_conf
    );

}