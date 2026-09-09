#include <vector>
#include <stdexcept>
#include <cassert>
#include <random>
#include <cstdint>
#include "cdl.hpp"

namespace channels::cdl {

    nrCDLChannel::nrCDLChannel() 
        : nrCDLChannel(Config{}) {
    }

    nrCDLChannel::nrCDLChannel(const Config& config) {

        configure(config);
    
    }

    const nrCDLChannel::Config& nrCDLChannel::config() const noexcept {

        return config_;
    }

    void nrCDLChannel::configure(const Config& config) {

        validateOrThrow(config);

        config_ = config;
        rng_.seed(config_.RandomStream.Seed);
        currentTime_ = config_.ChannelControl.InitialTime;

        initializeChannel();
    }

    void nrCDLChannel::initializeChannel() {

        // Step 1: configure antenna structure
        antenna::AntennaSystemConfig ant_sys_conf = 
            antenna::initializeAntennaStructure(
                config_.TransmitAntennaArraySetup,
                config_.ReceiveAntennaArraySetup
            );

        // Step 2: initialize delay profile
        // initializeDelayProfile();

        // Step 3: split LOS cluster and perform subclustering
        // buildClusters();

        // Step 4: generate initial phases
        // generateInitialPhases();

        // Step 5: compute ray coupling
        // coupleRays();

        // Step 6: initialize dual-mobility scatterer variables
        // initializeMobilityState();

        // Step 7: generate static CDL channel
        // generateChannelRealization();

        // Step 8: generate initial time-varying CDL channel

    }

    void nrCDLChannel::advance()
    {
        // Advance time
        //currentTime_ += timeStep_;

        // Evolve path gains according to:
        //
        // - Doppler spectrum
        // - fading distribution
        // - moving scatterers
        // - phase evolution
        // - whatever other time-dependent terms your CDL implementation uses

        //updatePathGains();
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
            delay.DelayProfile != pdp::PDP::CDL_D &&
            delay.DelayProfile != pdp::PDP::CDL_E) {
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