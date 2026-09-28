#include <vector>
#include <stdexcept>
#include <cassert>
#include "antenna_array/antenna.hpp"
#include "phase/initial_phase.hpp"
#include "cluster_profiles/profile.hpp"
#include "ray_coupling/coupling.hpp"
#include "mobility/mobility.hpp"
#include "response/generator.hpp"
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

        antenna::AntennaSystemConfig& AntennaSystem = config_.AntennaSystem;
        pdp::DelayProfileConfig& DelayProfile = config_.DelayProfile;
        ChannelControlConfig& ChannelControl = config_.ChannelControl;
        RandomStreamConfig& RandomStream = config_.RandomStream;
        mobility::MobilityConfig& Mobility = config_.Mobility;

        // Step 1: configure antenna structure
        antenna::geometry::initialize(AntennaSystem);

        // Step 2: initialize delay profile
        pdp::initializeDelayProfile(DelayProfile);

        // TODO: TBI
        // Step 3: split LOS cluster and perform subclustering
        // buildClusters(); 

        // Step 4: generate initial phases
        std::vector<float> phi = phase::generateInitialPhase(
            DelayProfile,
            RandomStream,
            ChannelControl,
            AntennaSystem
        );

        // Step 5: compute ray coupling
        std::vector<std::size_t> ray_coupling = ray::compute_coupling(
            DelayProfile,
            RandomStream
        );

        // Step 6: initialize dual-mobility scatterer variables
        mobility::compute_scatterer_variables(
            Mobility, 
            RandomStream, 
            DelayProfile
        );

        // Step 7: generate static CDL channel
        std::vector<float> Hstatic = response::generate_static_path_gains(
            DelayProfile, AntennaSystem, ray_coupling, phi);

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
            delay.ProfileName != pdp::PDP::CDL_D &&
            delay.ProfileName != pdp::PDP::CDL_E) {
            throw std::invalid_argument(
                "nrCDLChannel >>> KFactor is only valid for CDL_D and CDL_E."
            );
        }

        if (std::get<float>(mobility.MaximumDopplerShift) < 0.0f) {
            throw std::invalid_argument(
                "nrCDLChannel >>> MaximumDopplerShift must be nonnegative."
            );
        }

        if (mobility.Scatterer.MovingScattererProportion < 0.0f ||
            mobility.Scatterer.MovingScattererProportion > 1.0f) {
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