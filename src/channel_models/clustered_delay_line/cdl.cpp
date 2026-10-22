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

    Channel::Channel() 
        : Channel(Config{}) {
    }

    Channel::Channel(const Config& config) {

        configure(config);
    }

    const Channel::Config& Channel::config() const noexcept {

        return config_;
    }

    void Channel::configure(const Config& config) {

        validateOrThrow(config);

        config_ = config;
        rng_.seed(config_.RandomStream.Seed);
        currentTime_ = config_.ChannelControl.InitialTime;
        initialized_state_ = false;

        initial_phase_.clear();
        ray_coupling_.clear();

    }

    const std::vector<std::complex<float>>& Channel::operator()() {

        if (!initialized_state_) {

            initializeChannel();
            initialized_state_ = true;
        }
        else {

            advance();
        }

        // Generate and return the selected channel response
        const auto& output = response::generate_response_output(
            config_.PathGain,
            config_.OFDMResponse,
            config_.ChannelControl,
            config_.DelayProfile,
            config_.AntennaSystem,
            config_.Mobility
        );

        initialized_state_ = true;

        return output;
    }

    void Channel::initializeChannel() {

        antenna::AntennaSystemConfig& AntennaSystem = config_.AntennaSystem;
        pdp::DelayProfileConfig& DelayProfile = config_.DelayProfile;
        ChannelControlConfig& ChannelControl = config_.ChannelControl;
        RandomStreamConfig& RandomStream = config_.RandomStream;
        mobility::MobilityConfig& Mobility = config_.Mobility;
        response::PathGainConfig& PathGain = config_.PathGain;
        response::OFDMResponseConfig& OFDMResponse = config_.OFDMResponse;

        // Step 1: Configure antenna structure
        antenna::geometry::initialize(AntennaSystem);

        // Step 2: Initialize delay profile
        pdp::initializeDelayProfile(DelayProfile);

        // TODO: TBI
        // Step 3: split LOS cluster and perform subclustering
        // buildClusters(); 

        // Step 4: Generate initial phases
        initial_phase_ = phase::generateInitialPhase(
            DelayProfile,
            RandomStream,
            ChannelControl,
            AntennaSystem
        );

        // Step 5: Generate ray coupling
        ray_coupling_ = ray::compute_coupling(
            DelayProfile,
            RandomStream
        );

        // Step 6: Initialize mobility variables
        mobility::compute_scatterer_variables_los_nlos(
            Mobility, 
            RandomStream, 
            DelayProfile
        );

        // Step 7: Generate static path gains
        response::generate_static_path_gains(
            PathGain.LoSStaticGain, 
            PathGain.NLoSStaticGain, 
            ChannelControl, 
            DelayProfile, 
            AntennaSystem, 
            ray_coupling_, 
            initial_phase_
        );

        
    }

    void Channel::advance()
    {
        // Advance time
        // currentTime_ += timeStep_;

        // Evolve path gains according to:
        //
        // - Doppler spectrum
        // - fading distribution
        // - moving scatterers
        // - phase evolution
        // - whatever other time-dependent terms your CDL implementation uses

        //updatePathGains();

    }

    void Channel::validateOrThrow(const Config& config) {
        const auto& delay = config.DelayProfile;
        const auto& mobility = config.Mobility;
        const auto& control = config.ChannelControl;

        if (delay.DelaySpread <= 0.0f) {
            throw std::invalid_argument(
                "Channel >>> DelaySpread must be positive."
            );
        }

        if (delay.KFactor.has_value() &&
            delay.ProfileName != pdp::PDP::CDL_D &&
            delay.ProfileName != pdp::PDP::CDL_E) {
            throw std::invalid_argument(
                "Channel >>> KFactor is only valid for CDL_D and CDL_E."
            );
        }

        if (std::get<float>(mobility.MaximumDopplerShift) < 0.0f) {
            throw std::invalid_argument(
                "Channel >>> MaximumDopplerShift must be nonnegative."
            );
        }

        if (mobility.Scatterer.MovingScattererProportion < 0.0f ||
            mobility.Scatterer.MovingScattererProportion > 1.0f) {
            throw std::invalid_argument(
                "Channel >>> MovingScattererProportion must be between 0 and 1."
            );
        }

        if (control.SampleRate <= 0.0f) {
            throw std::invalid_argument(
                "Channel >>> SampleRate must be positive."
            );
        }

        if (control.SampleDensity <= 0.0f) {
            throw std::invalid_argument(
                "Channel >>> SampleDensity must be positive."
            );
        }

        if (control.NumTimeSamples == 0) {
            throw std::invalid_argument(
                "Channel >>> NumTimeSamples must be greater than zero."
            );
        }
    }

}