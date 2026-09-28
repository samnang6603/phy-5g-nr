#include "mobility.hpp"
#include "../cluster_profiles/profile.hpp"
#include <variant>

namespace channels::cdl::mobility {

    void compute_scatterer_variables(
        MobilityConfig& mobility_conf,
        RandomStreamConfig& randstream_conf,
        const pdp::DelayProfileConfig& pdp_conf
    ) {
        // Compute moving scatterer variables

        ScattererConfig& scatter_conf = mobility_conf.Scatterer;
        std::vector<bool>& states = scatter_conf.States;
        std::vector<float>& speed = scatter_conf.Speed;

        // TODO: get subcluster NLOS

        // Scatterer proportion dictates how many % of moving scatterers relative to
        // all clusters

        float p = scatter_conf.MovingScattererProportion;

        float v_scatt = scatter_conf.MaximumScattererSpeed;

        if (std::holds_alternative<float>(mobility_conf.MaximumDopplerShift)) {
            v_scatt = 0;
        }

        uint8_t M = NUMBER_OF_RAYS_PER_CLUSTER;
        uint8_t L = pdp_conf.ClusterTypes.size();
        std::size_t size = static_cast<std::size_t>(M)*L;

        auto& randomStream = randstream_conf.Stream;

        states.resize(size);
        speed.resize(size, 0.0f);

        // Use separate loops so each array consumes a contiguous
        // sequence of values from the random stream.
        for (auto v : states) {
            // Find out how many moving scatterers
            v = (static_cast<float>(randomStream.rand()) < p);
        }

        // For LOS cluster, there is no moving scatterers
        if (pdp_conf.HasLoS) {
            for (std::size_t i = 0; i < L; ++i) {
                states[i*M] = false;
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
}