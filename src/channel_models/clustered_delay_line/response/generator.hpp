#pragma once

#include <array>
#include "../antenna_array/antenna.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"

// Ray offset in 7.5 Table 7.5-3
constexpr std::array<float, NUMBER_OF_RAYS_PER_CLUSTER> RAY_OFFSET_ALPHA = {
    0.0447f,  -0.0447f,   0.1413f,  -0.1413f,   0.2492f,  
   -0.2492f,   0.3715f,  -0.3715f,   0.5129f,  -0.5129f,
    0.6797f,  -0.6797f,   0.8844f,  -0.8844f,  1.1481f,
   -1.1481f,   1.5195f,  -1.5195f,   2.1551f,  -2.1551f
};

namespace channels::cdl::response {

    std::vector<float> generate_static_path_gains(
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const std::vector<std::size_t>& coupling,
        const std::vector<float>& phi
    );
    //void generate_time_varying_path_gains();
    //void generate_ofdm_response();
    //void generate();
}