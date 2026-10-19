#pragma once
#include <array>
#include <variant>
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
        // Max doppler array can be 2 elem due to doppler from both rx and tx
        std::variant<float, std::array<float, 2>> MaximumDopplerShift{5.0f}; 
        ScattererConfig Scatterer;

        // User-Terminal direction of travel can be a 2 or 4-elements
        // if 2, then the values map to [Rx_azimuth, Rx_zenith]
        // if 4, then the values map to [Rx_azmith, Tx_azimuth, Rx_zenith, Tx_zenith]
        // azimuth is theta, zenith is phi
        std::vector<float> UTDirectionOfTravel = {0.0f, 90.0f};
    };

    void compute_scatterer_variables(
        MobilityConfig& mobility_conf,
        RandomStreamConfig& randstream_conf,
        const pdp::DelayProfileConfig& pdp_conf
    );


}