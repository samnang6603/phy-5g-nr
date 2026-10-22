#include <algorithm>
#include "mobility.hpp"
#include "../cluster_profiles/profile.hpp"
#include <variant>

namespace channels::cdl::mobility {

    void compute_scatterer_variables(
        MobilityConfig& mobility_conf,
        RandomStreamConfig& randstream_conf,
        const pdp::DelayProfileConfig& pdp_conf
    ) {
        // Compute moving scatterer variables.

        ScattererConfig& scatter_conf = mobility_conf.Scatterer;
        std::vector<uint8_t>& states = scatter_conf.States;
        std::vector<float>& speed = scatter_conf.Speed;

        // TODO: get subcluster NLOS

        // Scatterer proportion dictates how many % of moving scatterers relative to
        // all clusters.
        float p = scatter_conf.MovingScattererProportion;
        float v_scatt = scatter_conf.MaximumScattererSpeed;
        if (std::holds_alternative<float>(mobility_conf.MaximumDopplerShift)) {
            v_scatt = 0;
        }

        uint8_t M = NUM_RAYS_PER_CLUSTER;
        uint8_t L = pdp_conf.ClusterTypes.size();
        std::size_t size = static_cast<std::size_t>(M)*L;

        auto& randomStream = randstream_conf.Stream;

        states.resize(size);
        speed.resize(size, 0.0f);

        // Use separate loops so each array consumes a contiguous
        // sequence of values from the random stream.
        for (auto v : states) {
            // Find out how many moving scatterers
            // Bernoulli distribution of probability p. This specifies whether the 
            // scatters are moving (true) or not (false), for each ray in each cluster
            v = (static_cast<float>(randomStream.rand()) < p);
        }

        // For LOS cluster, there is no moving scatterers
        if (pdp_conf.HasLoS) {
            for (std::size_t i = 0; i < L; ++i) {
                states[i*M] = 0;
            } 
        }

        // Generate RV D (size [N M]) between -v_scatt to v_scatt
        // randD is (0,1) so (0,1)*2-1 = (-1,1). (-1,1)*v_scatt = (-v_scatt,v_scatt)
        if (v_scatt) {
            for (auto& v : speed) {
                v = (static_cast<float>(randomStream.rand())*2 - 1)*v_scatt;
            }
        }

        // To be implemented later, available subcluster scenario

    }

    void compute_scatterer_variables_los_nlos(
        MobilityConfig& mobility_conf,
        RandomStreamConfig& randstream_conf,
        const pdp::DelayProfileConfig& pdp_conf
    ) {
        // Compute moving scatterer variables

        ScattererConfig& scatter_conf = mobility_conf.Scatterer;
        std::vector<uint8_t>& los_states = scatter_conf.LoSStates;
        std::vector<float>& los_speed = scatter_conf.LoSSpeed;
        std::vector<uint8_t>& nlos_states = scatter_conf.NLoSStates;
        std::vector<float>& nlos_speed = scatter_conf.NLoSSpeed;

        // TODO: get subcluster NLOS

        // Scatterer proportion dictates how many % of moving scatterers relative to
        // all clusters.
        float p = scatter_conf.MovingScattererProportion;
        float v_scatt = scatter_conf.MaximumScattererSpeed;
        if (std::holds_alternative<float>(mobility_conf.MaximumDopplerShift)) {
            v_scatt = 0;
        }

        uint8_t M = NUM_RAYS_PER_CLUSTER;
        uint8_t L = pdp_conf.ClusterTypes.size();
        uint8_t L1 = L - static_cast<uint8_t>(pdp_conf.HasLoS);
        std::size_t nlos_size = static_cast<std::size_t>(M)*(L1);
        nlos_states.resize(nlos_size);
        nlos_speed.resize(nlos_size, 0.0f);
        los_states.resize(M);
        los_speed.resize(M);

        auto& randomStream = randstream_conf.Stream;

        uint8_t* stp = nlos_states.data();
        float* svp = nlos_speed.data();

        // This implementation is weird because it's (for now) intended to match MATLAB golden model
        // Channel() bit-exact. As soon as work is done and all cases are corrected, a simpler implementation
        // (no skipping voiding output just to consume random stream) can be applied.
        if (pdp_conf.HasLoS) {

            std::fill(los_speed.begin(), los_speed.end() ,0.0f);
            std::fill(los_states.begin(), los_states.end(), static_cast<uint8_t>(0));

            for (std::size_t m = 0; m < M; ++m) {

                // Consume the first sample of this column
                (void)randomStream.rand();

                // Store the remaining M-1 samples
                for (std::size_t l = 0; l < L1; ++l) {

                    *stp++ = static_cast<uint8_t>(randomStream.rand() < p);
                }
            } 

            // Generate RV D (size [N M]) between -v_scatt to v_scatt
            // randD is (0,1) so (0,1)*2-1 = (-1,1). (-1,1)*v_scatt = (-v_scatt,v_scatt)
            if (v_scatt) {
                for (std::size_t m = 0; m < M; ++m) {

                // Consume the first sample of this column
                (void)randomStream.rand();

                    // Store the remaining M-1 samples
                    for (std::size_t l = 0; l < L1; ++l) {
                        *svp++ = static_cast<float>(randomStream.rand()*2 - 1)*v_scatt;
                    }
                } 
            }

            // To be implemented later, available subcluster scenario

        } else {

            for (std::size_t m = 0; m < M; ++m) {
                for (std::size_t l = 0; l < L; ++l) {
                    *stp++ = static_cast<uint8_t>(randomStream.rand() < p);
                }
            } 

            // Generate RV D (size [N M]) between -v_scatt to v_scatt
            // randD is (0,1) so (0,1)*2-1 = (-1,1). (-1,1)*v_scatt = (-v_scatt,v_scatt)
            if (v_scatt) {
                for (std::size_t m = 0; m < M; ++m) {

                    for (std::size_t l = 0; l < L; ++l) {
                        *svp++ = static_cast<float>(randomStream.rand()*2 - 1)*v_scatt;
                    }
                } 
            }
        }
    }
}