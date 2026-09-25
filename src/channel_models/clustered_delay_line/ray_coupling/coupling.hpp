#pragma once

#include <vector>
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"

namespace channels::cdl::ray {

    std::vector<std::size_t> compute_coupling(
        pdp::DelayProfileConfig& pdp_conf,
        RandomStreamConfig& rstream_conf
    );
    
}