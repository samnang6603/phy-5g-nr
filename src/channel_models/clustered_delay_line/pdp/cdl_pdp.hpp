#pragma once

#include <optional>
#include "cdl_pdp_lut.hpp"

namespace channels::cdl::pdp {

    enum class PDP { CDL_A, CDL_B, CDL_C, CDL_D, CDL_E };

    struct MeanAnglesList {
        float AoD = 0.0f;
        float AoA = 0.0f;
        float ZoD = 0.0f;
        float ZoA = 0.0f;
    };

    struct DelayProfileConfig {
        PDP DelayProfile = PDP::CDL_A;
        float DelaySpread = 3E-8f;

        // nullopt means disabled
        std::optional<float> KFactor;

        // nullopt means angle scaling disabled
        std::optional<MeanAnglesList> MeanAngles;
    };



}

