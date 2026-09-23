#pragma once

#include <vector>
#include "cdl_pdp_lut.hpp"
#include <optional>

namespace channels::cdl::pdp {

    enum class PDP { CDL_A, CDL_B, CDL_C, CDL_D, CDL_E };

    struct MeanAnglesList {
        float AoD = 0.0f;
        float AoA = 0.0f;
        float ZoD = 0.0f;
        float ZoA = 0.0f;
    };

    struct AngleSpreadsList {
        float C_ASD;
        float C_ASA;
        float C_ZSD;
        float C_ZSA;
    };

    struct DelayProfileConfig {
        PDP ProfileName = PDP::CDL_A;
        bool HasLoS = false;
        float DelaySpread = 3E-8f;
        std::vector<CDLCluster> Table;

        // nullopt means disabled
        std::optional<float> KFactor;

        // nullopt means angle scaling disabled
        std::optional<MeanAnglesList> MeanAngles;

        AngleSpreadsList AngleSpreads;
        float XPR;

    };

    void initializeDelayProfile(DelayProfileConfig& pdp_conf);



}

