#pragma once

#include <vector>
#include "../cluster_profiles/profile.hpp"
#include "../../channel_models_common.hpp"

namespace channels::cdl::mobility {

    struct ScattererConfig {
        float MaximumScattererSpeed = 5.0f;
        float MovingScattererProportion = 0.2f;
        std::vector<float> Speed;
        std::vector<bool> States;
    };

    struct MobilityConfig {
        float MaximumDopplerShift = 5.0f;
        ScattererConfig Scatterer;
    };

    void compute_scatterer_variables(
        MobilityConfig& mobility_conf,
        RandomStreamConfig& randstream_conf,
        const pdp::DelayProfileConfig& pdp_conf
    );


}