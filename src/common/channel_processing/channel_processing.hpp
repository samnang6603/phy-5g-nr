#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>
#include <cmath>
#include <complex>

namespace common::channel_processing {

    inline int estimate_delay(
        const std::vector<std::complex<float>>& path_gains,
        const std::vector<float>& path_filters,
        const std::size_t num_time_samp,
        const std::size_t num_paths,
        const std::size_t num_Tx,
        const std::size_t num_Rx
    ) {
        // Perfect delay estimation of channel by finding the peak of the channel impulse
        // response. The impulse response is averaged across all channel time samples (snapshots)
        // before being summed across all transmit and then receive antennas.
        // Assumes path_gains has the flat arrangement of Ncs x Npath x NTx x NRx.
        //          
        //                 ___ | ___  ___               |
        //  delay = argmax \   | \    \   g_avg_p*f_p[k]|
        //            k    /__ | /__  /__               |
        //                  Rx |  Tx  paths             |
        //

        const std::size_t num_path_gains = path_gains.size();

        // Average the path gains array across all time elements (channel snapshot)
        std::vector<std::complex<float>> avg_gain_cs(num_path_gains/num_time_samp); 
        std::complex<float>* cs = avg_gain_cs.data();
        for (std::size_t i = 0; i < num_path_gains; i += num_time_samp) {

            for (std::size_t k = 0; k < num_time_samp; ++k) {
                *cs += path_gains[i + k]; 
            }
            *(cs++) /= num_time_samp;
        }

        const std::size_t filt_len_per_path = path_filters.size()/num_paths;
        std::vector<std::complex<float>> h(num_Rx*num_Tx*filt_len_per_path);
        //std::complex<float>* h_ptr = h.data();

        // TODO: Reorganize path_filters in channel_filter.cpp so that path
        //       is the contiguous dimension. This will improve memory access
        //       in the innermost loop and enable more efficient SIMD vectorization.
        //
        // TODO: After validating the optimized implementation across multiple
        //       simulations and edge cases, remove the debug variables: 
        //       tx_mag_sum, stx, h_ptr, and h.
        //
        // Compute the channel response magnitude summed across all Tx/Rx antenna pairs.
        // For each filter tap (k):
        //   1. Sum the complex contributions from all paths (p) for each Tx/Rx pair.
        //   2. Coherently sum the complex responses across Tx antennas (s).
        //   3. Take the magnitude of the combined Tx response for each Rx antenna (u).
        //   4. Sum these magnitudes across Rx antennas to obtain the final tap value.
        //
        // Note: Tx responses are summed before taking the magnitude (coherent sum),
        //       while Rx magnitudes are summed afterward (non-coherent sum).
        std::complex<float> gain_filter_accum;
        std::vector<float> tx_mag_sum(filt_len_per_path*num_Rx);
        float* stx = tx_mag_sum.data();
        std::complex<float> stx_accum;
        std::vector<float> tx_rx_mag_sum(filt_len_per_path);
        float* stxrx = tx_rx_mag_sum.data();
        float stxrx_accum;
        float tmp;

        for (std::size_t k = 0; k < filt_len_per_path; ++k) {

            // Accumulate the magnitudes across all Rx antennas for tap k.
            stxrx_accum = 0;
            for (std::size_t u = 0; u < num_Rx; ++u) {

                // Coherently accumulate the Tx responses for Rx antenna u.
                stx_accum = 0;
                const std::size_t u_stride = u*num_Tx*num_paths;
                for (std::size_t s = 0; s < num_Tx; ++s) {

                    // Sum the complex path contributions for Tx/Rx pair.
                    gain_filter_accum = 0;
                    const std::size_t us_stride = u_stride + s*num_paths;
                    for (std::size_t p = 0; p < num_paths; ++p) {
                        const std::size_t pf_idx = k + filt_len_per_path*p;
                        gain_filter_accum += path_filters[pf_idx]*avg_gain_cs[us_stride + p];
                    }

                    //*h_ptr++ = gain_filter_accum;
                    stx_accum += gain_filter_accum;
                }

                // Convert the coherent Tx sum to magnitude, then accumulate
                // the magnitudes across Rx antennas.
                tmp = std::abs(stx_accum);
                *stx++ = tmp;
                stxrx_accum += tmp;
            }

            // Store the combined Tx/Rx magnitude for filter tap k.
            *stxrx++ = stxrx_accum;
        }

        auto peak_idx_iter = std::max_element(tx_rx_mag_sum.begin(), tx_rx_mag_sum.end());
        const int delay = peak_idx_iter - tx_rx_mag_sum.begin();
        return delay;
    }

} // namespace common::channel_processing