#pragma once

#include <random>
#include "cluster_profiles/profile.hpp"
#include "antenna_array/antenna.hpp"
#include "mobility/mobility.hpp"
#include "response/generator.hpp"
#include "../channel_models_common.hpp"

namespace channels::cdl {

    class Channel {
    public:

        struct Config {
            pdp::DelayProfileConfig DelayProfile;
            mobility::MobilityConfig Mobility;
            antenna::AntennaSystemConfig AntennaSystem;
            RandomStreamConfig RandomStream;
            ChannelControlConfig ChannelControl;
            response::PathGainConfig PathGain;
            response::OFDMResponseConfig OFDMResponse;
        };

        Channel();

        explicit Channel(const Config& config);

        const Config& config() const noexcept;
        
        void configure(const Config& config);

        const std::vector<std::complex<float>>& operator()();

        void advance();

    private:

        void initializeChannel();

        static void validateOrThrow(const Config& config);

    private:

        Config config_;

        // Random number generator
        std::mt19937 rng_;

        // Simulation time
        float currentTime_ = 0.0f;

        // Channel initialization state
        bool initialized_state_ = false;

        // Persistent channel realization variables
        std::vector<float> initial_phase_;
        std::vector<std::size_t> ray_coupling_;

        // Reference to selected response vector
        const std::vector<std::complex<float>> output_;
    };
}