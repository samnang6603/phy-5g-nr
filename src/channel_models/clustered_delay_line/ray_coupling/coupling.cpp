#include <algorithm>
#include <numeric>
#include <vector>
#include "coupling.hpp"
#include "../../channel_models_common.hpp"

// Coupling permutation
// 1. AoD -> AoA coupling
// 2: AoD -> ZoA coupling
// 3. AoD -> ZoD coupling
constexpr uint16_t COUPLING_PERMUTATIONS = 3;

namespace channels::cdl::ray {

    static void adjust_NLOS_coupling(
        std::vector<std::size_t>& coupling,
        const std::vector<PropagationCondition>& cluster_types,
        const bool hasLoS
    );

    static std::vector<std::size_t> sort_ray_indices(
        const std::vector<float>& x,
        const std::size_t L,
        const std::size_t M,
        const std::size_t P
    );

    std::vector<std::size_t> compute_coupling(
        pdp::DelayProfileConfig& pdp_conf,
        RandomStreamConfig& rstream_conf
    ) {
        // Calculate ray coupling according to TR 38.901 section 7.5 step 8

        std::vector<PropagationCondition>& cluster_types = pdp_conf.ClusterTypes;
        const std::uint8_t M = NUMBER_OF_RAYS_PER_CLUSTER;
        const std::uint8_t L = cluster_types.size();

        // First L*M plane is AoD -> AoA coupling
        // Second L*M plane is AoD -> ZoA coupling
        // Third L*M plane is AoD -> ZoD coupling
        const std::uint16_t coupling_size = L*M*COUPLING_PERMUTATIONS;
        std::vector<float> tmp(coupling_size);

        auto& randomStream = rstream_conf.Stream;

        for (auto& v : tmp) {
            v = static_cast<float>(randomStream.rand());
        }

        // Sort weakest to strongest coupling within ray per cluster per angle
        std::vector<std::size_t> indices = sort_ray_indices(tmp, L, M, COUPLING_PERMUTATIONS);

        /*
         Rearranging ray coupling to preserve directional consistency:
         -------------------------------------------------------------
         While ray angles (AoD, AoA, ZoD, ZoA) are randomly generated, 
         we want to avoid complete chaos. 
         Specifically, we preserve the sorted ZoD ↔ ZoA relationship across 
         rays for each cluster.
        
         This maintains a logical pairing between departure and arrival 
         directions (e.g., rays leaving high tend to arrive high), which 
         reflects realistic scattering behavior.
        
         Without this, the model would assign arrival angles that have no 
         directional correlation with their departure, breaking spatial 
         symmetry, phase continuity, and any realistic MIMO behavior.
        
         NOTE: Even if RayCoupling = 'Random', this logic still maintains 
         partial directional consistency by sorting ZoD -> ZoA mappings.
        */

        std::vector<std::size_t> coupling(coupling_size,0);

        // Copy AoD -> AoA
        std::copy(indices.data(),
          indices.data() + L*M,
          coupling.data()
        );

        const std::size_t AoD2ZoD_start_idx = (COUPLING_PERMUTATIONS - 1)*L*M;
        const std::size_t AoD2ZoA_start_idx = L*M;

        // Copy AoD -> ZoD
        std::copy(indices.data() + AoD2ZoD_start_idx,
          indices.data() + COUPLING_PERMUTATIONS*L*M,
          coupling.data() + AoD2ZoD_start_idx
        );

        std::vector<std::size_t> cluster_sorted_indices(M);
        std::vector<std::size_t> this_cluster(M);
        for (std::size_t l = 0; l < L; ++l) {
        
            std::iota(cluster_sorted_indices.begin(), cluster_sorted_indices.end(), 0);

            const auto* zod = indices.data() + AoD2ZoD_start_idx + l;

            // sort the AoD->ZoD in ascending order and get only the index
            std::sort(cluster_sorted_indices.begin(), cluster_sorted_indices.end(),
                [zod, L](std::size_t a, std::size_t b) {
                return zod[a*L] < zod[b*L]; 
            });

            // map the sorted index to rearrange AoD->ZoA according to AoD->ZoD
            auto* dst = coupling.data() + AoD2ZoA_start_idx + l;
            auto* zoa = indices.data() + AoD2ZoA_start_idx + l;
            for (std::size_t m = 0; m < M; ++m) {
                //coupling[AoD2ZoA_start_idx + l + L*m] = indices[AoD2ZoA_start_idx + l + L*cluster_sorted_indices[m]];
                *dst = *(zoa + L*cluster_sorted_indices[m]);
                dst += L;
            }
        }

        // TODO: Subcluster to be implemented later here

        adjust_NLOS_coupling(coupling, cluster_types, pdp_conf.HasLoS);

        return coupling;

    }

    static void adjust_NLOS_coupling(
        std::vector<std::size_t>& coupling,
        const std::vector<PropagationCondition>& cluster_types,
        const bool hasLoS
    ) {

        /*
        Each cluster has local indices [1 to M] for rays.
        We now remap those to global indices across all clusters of this type
        Think of each column as a ray index inside its own cluster.
        To globalize it, we stride by cluster. That means:
        For each ray index i in cluster k:
        global_index = [i - 1]*stride + k
        */

        /* TODO: To be implemented later when there're subclusters 
        std::vector<std::size_t> ind(cluster_types.size() - 1); 
        std::size_t stride = 0;
        for (auto& v : cluster_types) {
            if (v == PropagationCondition::NLOS) {
                ++stride;
            }
        } */

        /* TODO: To be deleted when subcluster is implemented */
        const std::size_t L = cluster_types.size();
        const std::size_t stride = hasLoS ? L - 1 : L;
        const std::size_t start = static_cast<std::size_t>(hasLoS);
        for (std::size_t l = start; l < L; ++l) {
            for (std::size_t m = 0; m < NUMBER_OF_RAYS_PER_CLUSTER; ++m) {
                for (std::size_t a = 0; a < COUPLING_PERMUTATIONS; ++a) {
                    const std::size_t this_idx = l + L*m + L*NUMBER_OF_RAYS_PER_CLUSTER*a;
                    auto& v = coupling[this_idx];
                    v = v*stride + l;
                }
            }
        }
    }

    static std::vector<std::size_t> sort_ray_indices(
        const std::vector<float>& x,
        const std::size_t L,
        const std::size_t M,
        const std::size_t P
    ) {
        // Sort ray indices on x and return the flat buffer of indices

        std::vector<std::size_t> ind_sorted(L*M*P);

        // Scratch index buffer
        std::vector<std::size_t> idx(M);

        for (std::size_t p = 0; p < P; ++p) {

            const std::size_t this_angle = p*L*M;

            for (std::size_t l = 0; l < L; ++l) {

                std::iota(idx.begin(), idx.end(),0);

                const float* base = x.data() + l + this_angle;

                std::sort(idx.begin(), idx.end(),
                [base, L] (std::size_t a, std::size_t b) {
                    return base[a*L] < base[b*L];
                });

                // Store using column-major layout
                for (std::size_t m = 0; m < M; ++m) {
                    ind_sorted[this_angle + l + m*L] = idx[m];
                }
            }
        }

        return ind_sorted;
    }

} // namespace channels::cdl::ray


