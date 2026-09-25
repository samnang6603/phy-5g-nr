#pragma once

#include <random>
#include "cluster_profiles/profile.hpp"
#include "antenna_array/antenna.hpp"
#include "mobility/mobility.hpp"
#include "../channel_models_common.hpp"

namespace channels::cdl {

    class nrCDLChannel {
    public:

        struct Config {
            pdp::DelayProfileConfig DelayProfile;
            mobility::MobilityConfig Mobility;
            antenna::AntennaSystemConfig AntennaSystem;
            RandomStreamConfig RandomStream;
            ChannelControlConfig ChannelControl;
        };

        nrCDLChannel();

        explicit nrCDLChannel(const Config& config);

        const Config& config() const noexcept;
        
        void configure(const Config& config);

        void advance();

    private:
        void initializeChannel();
        static void validateOrThrow(const Config& config);

    private:
        Config config_;
        std::mt19937 rng_;
        float currentTime_ = 0.0f;

    };
}