#pragma once

#include <cstdint>
#include <optional>
#include <random>
#include "pdp/cdl_pdp.hpp"
#include "antenna_structure/antenna_structure.hpp"
#include "../matlab_mt19937.hpp"

constexpr std::uint8_t NUMBER_OF_RAYS = 20;
constexpr std::uint64_t SPEED_OF_LIGHT = 299792458;

namespace channels::cdl {

    enum class ChannelResponseOutputType {
        PATH_GAINS,
        OFDM_RESPONSE
    };

    class nrCDLChannel {
    public:

        struct MobilityConfig {
            float MaximumDopplerShift = 5.0f;
            float MovingScattererProportion = 0.2f;
        };

        struct RandomStreamConfig {
            std::uint32_t Seed = 73;
            //MATLABMT19937 Stream;
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

        struct Config {
            pdp::DelayProfileConfig DelayProfile;
            MobilityConfig Mobility;
            RandomStreamConfig RandomStream;
            ChannelControlConfig ChannelControl;
            antenna::AntennaArrayConfig TransmitAntennaArraySetup;
            antenna::AntennaArrayConfig ReceiveAntennaArraySetup;
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