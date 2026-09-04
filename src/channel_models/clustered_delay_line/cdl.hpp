#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../nr5g_common.hpp"
#include "antenna_structure/antenna_structure.hpp"

static constexpr std::uint8_t NUMBER_OF_RAYS = 20;

namespace channels::cdl {

    enum class PDP { CDL_A, CDL_B, CDL_C, CDL_D, CDL_E };

    enum class ChannelResponseOutputType {
        PATH_GAINS,
        OFDM_RESPONSE
    };

    struct MeanAnglesList {
        float AoD = 0.0f;
        float AoA = 0.0f;
        float ZoD = 0.0f;
        float ZoA = 0.0f;
    };

    class nrCDLChannel {
    public:
        struct DelayProfileConfig {
            PDP DelayProfile = PDP::CDL_A;
            float DelaySpread = 3E-8f;

            // nullopt means disabled
            std::optional<float> KFactor;

            // nullopt means angle scaling disabled
            std::optional<MeanAnglesList> MeanAngles;
        };

        struct MobilityConfig {
            float MaximumDopplerShift = 5.0f;
            float MovingScattererProportion = 0.2f;
        };

        struct RandomStreamConfig {
            std::uint32_t Seed = 1;
        };

        struct ChannelFilteringConfig {
            uint64_t FilterDelay = 7;               // samples
            float StopbandAttenuation = 70.0f;     // dB
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
        };

        struct Config {
            DelayProfileConfig DelayProfile;
            MobilityConfig Mobility;
            RandomStreamConfig RandomStream;
            ChannelControlConfig ChannelControl;
        };

        nrCDLChannel();

        explicit nrCDLChannel(const Config& config);

        const Config& config() const noexcept;
        
        void configure(const Config& config);

    private:
        Config config_;
        std::mt19937 rng_;
        float currentTime_ = 0.0f;

    private:
        static void validateOrThrow(const Config& config);
    };
}