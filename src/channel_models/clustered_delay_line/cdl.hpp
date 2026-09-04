#pragma once

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <optional>

#include "../../nr5g_common.hpp"
#include "antenna_structure/antenna_structure.hpp"

/******************* Constants ***************************/
static constexpr uint8_t NUMBER_OF_RAYS = 20;

/***************** Type Definitions **********************/
namespace channels::cdl {
    
    enum class PDP { CDL_A, CDL_B, CDL_C, CDL_D, CDL_E };

    struct MeanAnglesList {
        float AoD = 0.0f; // mean azimuth of departure angles after scaling
        float AoA = 0.0f; // mean azimuth of arrival angles after scaling
        float ZoD = 0.0f; // mean zenith of departure angles after scaling
        float ZoA = 0.0f; // mean zenith of arrival angles after scaling
    };

    class nrCDLChannel {
    public:
        struct KFactorScalingConfig {
            float KFactor = 9.0f;
        };

        struct AngleScalingConfig {
            MeanAnglesList MeanAngles{};
        };

        struct Config {
            PDP DelayProfile = PDP::CDL_A;
            float DelaySpread = 3E-8f;

            // If nullopt, K-factor scaling is disabled
            std::optional<KFactorScalingConfig> KFactorScaling;

            // If nullopt, angle scaling is disabled
            std::optional<AngleScalingConfig> AngleScaling;
        };

    public:
        nrCDLChannel()
            : config_{} {
            validateOrThrow(config_);
        }

        explicit nrCDLChannel(Config config)
            : config_{config} {
            validateOrThrow(config_);
        }

        const Config& config() const noexcept {
            return config_;
        }

        void configure(Config config) {
            validateOrThrow(config);
            config_ = config;
        }

    private:
        Config config_;

    private:
        static constexpr bool supportsKFactorScaling(PDP pdp) noexcept {
            return pdp == PDP::CDL_D || pdp == PDP::CDL_E;
        }

        static void validateOrThrow(const Config& config) {
            if (config.DelaySpread <= 0.0f) {
                throw std::invalid_argument("nrCDLChannel >>> DelaySpread must be positive.");
            }

            if (config.KFactorScaling.has_value() &&
                !supportsKFactorScaling(config.DelayProfile)) {
                throw std::invalid_argument(
                    "nrCDLChannel >>> KFactorScaling is only valid for CDL_D and CDL_E."
                );
            }
        }
    };

}