#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include "matlab_mt19937.hpp"

constexpr std::uint8_t NUMBER_OF_RAYS_PER_CLUSTER = 20;
constexpr float SPEED_OF_LIGHT = 299792458.0f;

namespace channels {

    enum class ChannelResponseOutputType {
        PATH_GAINS,
        OFDM_RESPONSE
    };

    enum class PropagationCondition {
        LOS,
        NLOS,
        SUBCLUSTER_NLOS,
        NUMBER_OF_PROPAGATION_CONDITIONS 
    };

    struct ChannelFilteringConfig {
        uint64_t FilterDelay = 7;               // samples
        float StopbandAttenuation = 70.0f;      // dB
        float MaxFractionalDelayError = 0.01f;
    };

    struct ChannelControlConfig {
        float SampleRate = 30720000.0f;
        float InitialTime = 0.0f;
        float SampleDensity = 64.0f;

        bool NormalizeChannelOutput = true;
        bool NormalizePathGains = true;

        ChannelResponseOutputType ChannelResponseOutput =
        ChannelResponseOutputType::PATH_GAINS;

        std::uint64_t NumTimeSamples = 30720;

        // nullopt means channel filtering disabled
        std::optional<ChannelFilteringConfig> ChannelFiltering;

        std::string InitialPhase = "Random";
        std::string RayCoupling  = "Random";

    };

    struct RandomStreamConfig {
        std::uint32_t Seed;
        MATLABMT19937 Stream;

        RandomStreamConfig()
            : Seed(73),
            Stream(Seed)
        {}
    };
}