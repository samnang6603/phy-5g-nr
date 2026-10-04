#pragma once

#include <vector>
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"

/***************** Type Definitions *****************/
namespace channels::cdl::ray {

// Coupling permutation
// 1. AoD -> AoA coupling
// 2: AoD -> ZoA coupling
// 3. AoD -> ZoD coupling
constexpr uint16_t COUPLING_PERMUTATIONS = 3;
constexpr uint16_t AOD_TO_AOA_COUPLING_PLANE = 1;
constexpr uint16_t AOD_TO_ZOA_COUPLING_PLANE = 2;
constexpr uint16_t AOD_TO_ZOD_COUPLING_PLANE = 3;

}

/*************** Function Declarations ****************/
namespace channels::cdl::ray {

    std::vector<std::size_t> compute_coupling(
        pdp::DelayProfileConfig& pdp_conf,
        RandomStreamConfig& rstream_conf
    );
    
}