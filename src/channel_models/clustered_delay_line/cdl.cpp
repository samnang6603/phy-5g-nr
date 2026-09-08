#include <vector>
#include <stdexcept>
#include <cassert>
#include <random>
#include <cstdint>
#include "cdl.hpp"

namespace channels::cdl {

    nrCDLChannel::nrCDLChannel() : nrCDLChannel(Config{}) {
    }

    nrCDLChannel::nrCDLChannel(const Config& config)
        : config_{config},
          rng_{config.RandomStream.Seed},
          currentTime_{config.ChannelControl.InitialTime} {

        validateOrThrow(config_);
    }

    const nrCDLChannel::Config& nrCDLChannel::config() const noexcept {

        return config_;
    }

    void nrCDLChannel::configure(Config& config) {

        validateOrThrow(config);

        config_ = config;
        rng_.seed(config_.RandomStream.Seed);
        currentTime_ = config_.ChannelControl.InitialTime;

        antenna::arrangeStructure(config.AntennaArraySetup);
    }

    void nrCDLChannel::validateOrThrow(const Config& config) {
        const auto& delay = config.DelayProfile;
        const auto& mobility = config.Mobility;
        const auto& control = config.ChannelControl;

        if (delay.DelaySpread <= 0.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> DelaySpread must be positive."
            );
        }

        if (delay.KFactor.has_value() &&
            delay.DelayProfile != PDP::CDL_D &&
            delay.DelayProfile != PDP::CDL_E) {
            throw std::invalid_argument(
                "nrCDLChannel >>> KFactor is only valid for CDL_D and CDL_E."
            );
        }

        if (mobility.MaximumDopplerShift < 0.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> MaximumDopplerShift must be nonnegative."
            );
        }

        if (mobility.MovingScattererProportion < 0.0f ||
            mobility.MovingScattererProportion > 1.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> MovingScattererProportion must be between 0 and 1."
            );
        }

        if (control.SampleRate <= 0.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> SampleRate must be positive."
            );
        }

        if (control.SampleDensity <= 0.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> SampleDensity must be positive."
            );
        }

        if (control.NumTimeSamples == 0) {
            throw std::invalid_argument(
                "nrCDLChannel >>> NumTimeSamples must be greater than zero."
            );
        }
    }

}