#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include "channel_filter.hpp"
#include "../../utils/utils.hpp"

namespace channels::filter {

    static float compute_kaiser_beta(const float atten);
    static std::array<float, FILTER_PROTOTYPE_LENGTH> design_fir_kaiser_interpolator();

    using Phase = std::array<float, FILTER_PROTOTYPE_LENGTH>;
    using PhaseBank = std::array<Phase, INTERPOLATION_FACTOR>;

    // MATLAB equivalent: E = reshape(polyFIR, 50, 16)
    // E(p+1,k+1) == polyFIR(p + 50*k + 1), using MATLAB 1-based indexing.
    PhaseBank polyphase(const std::array<float, FILTER_PROTOTYPE_LENGTH>& b) {
        PhaseBank phases{};
        for (std::size_t p = 0; p < INTERPOLATION_FACTOR; ++p) {
            for (std::size_t k = 0; k < FILTER_LENGTH; ++k) {
                phases[p][k] = b[p + k*INTERPOLATION_FACTOR];
            }
        }
        return phases;
    }

    void design_channel_filter(
        ChannelFilterConfig& chan_filt_conf,
        const std::vector<float> path_delays
    ) {
        // Design channel polyphase fractional interpolation filter

        const std::size_t num_paths = path_delays.size();
        const float fs = chan_filt_conf.fs;

        //std::vector<float> frac_delay(num_paths);
        std::vector<std::size_t> phase_idx(num_paths);
        std::vector<float> frac_delay_range(INTERPOLATION_FACTOR + 1);
        std::vector<int> tap_range(num_paths*2);
        std::vector<std::size_t> coeff_range(num_paths);
        
        for (std::size_t i = 0; i < num_paths; ++i) {
            const float pf = path_delays[i]*fs;
            const float frac_delay = pf - std::floor(pf);

            for (std::size_t k = 0; k <= INTERPOLATION_FACTOR; ++k) {
                frac_delay_range[k] = std::fabs(frac_delay - INTERPOLATION_RANGE[k]);
            }

            phase_idx[i] = std::min_element(
                frac_delay_range.begin(), 
                frac_delay_range.end()) - frac_delay_range.begin();

            const bool is_not_endpt = !((phase_idx[i] == 0) || (phase_idx[i] == INTERPOLATION_FACTOR));
            const int is_not_endpt_int = static_cast<int>(is_not_endpt);
            const int quantized_last = static_cast<int>(phase_idx[i] == INTERPOLATION_FACTOR);

            // Calculate integer delay part
            const int int_delay = std::floor(pf) + quantized_last;

            // Tap index range, excluding the path with fracDelay closest to 0 or 1
            // - Tap delay refer to propagation environment (CDL/TDL). It affects actual
            // multipath spacing. Each tap is one cluster/echo path with a certain delay
            // comeing from PDP. Unit is usually nanosecond.
            // - Channel filter delay refer to FIR filter (interpolation/smoothing). It
            // affects only the filtering mechanics and output alignment. Unit is
            // usually in samples.
            // The + FILTER_DELAY is updating integer delays
            tap_range[i] = int_delay + is_not_endpt_int*(-FILTER_LENGTH_HALF + 1) + FILTER_DELAY;
            coeff_range[i] = tap_range[i] + is_not_endpt_int*(FILTER_LENGTH - 1) + 1;
            tap_range[i + num_paths] = int_delay + is_not_endpt_int*FILTER_LENGTH_HALF - 1;
        }
        
        const std::size_t coeff_max_size = *std::max_element(coeff_range.begin(), coeff_range.end());
        std::vector<float>& h_interp = chan_filt_conf.path_filter_coeffs;
        h_interp.resize(num_paths*coeff_max_size);
        const std::array<float, FILTER_PROTOTYPE_LENGTH> polyphase_fir = design_fir_kaiser_interpolator();
        const float* pfir;

        for (std::size_t i = 0; i < num_paths; ++i) {

            const std::size_t coeff_start_idx = tap_range[i];
            const bool is_endpt = (phase_idx[i] == 0) || (phase_idx[i] == INTERPOLATION_FACTOR);

            if (is_endpt) {
                // don't use this if it's at end point
                //h_chan_filt[i*FILTER_LENGTH + FILTER_LENGTH_HALF - 1] = 1.0f;
                h_interp[i*coeff_max_size + coeff_start_idx] = 1.0f;
                continue;
            }

            pfir = polyphase_fir.data() + phase_idx[i];
            
            //for (std::size_t k = 0; k < FILTER_LENGTH; ++k) {
            //    h_chan_filt[i*FILTER_LENGTH + k] = *(pfir + k*INTERPOLATION_FACTOR);
            //}

            for (std::size_t k = 0; k < FILTER_LENGTH; ++k) {
                h_interp[i*coeff_max_size + coeff_start_idx + k] = *(pfir + k*INTERPOLATION_FACTOR);
            }
        }
    }

    static std::array<float, FILTER_PROTOTYPE_LENGTH> design_fir_kaiser_interpolator() {

        constexpr uint16_t L = INTERPOLATION_FACTOR;
        constexpr uint16_t M = 1; // redundant but to stress that this is a pure interpolator
        constexpr float R = static_cast<float>(L)/M;
        constexpr float R_inv = 1.0f/R;
        constexpr float center = static_cast<float>(IMPULSE_RESPONSE_CENTER);
        constexpr float center_inv = 1.0f/center;

        const float beta = compute_kaiser_beta(STOPBAND_ATTENUATION);

        const float I0_beta = utils::float_arithmetic::besselI0(beta);
        const float I0_beta_inv = 1.0f/I0_beta;

        std::array<float, FILTER_PROTOTYPE_LENGTH> h{}; // initialize Nyquist zeros to 0

        // Process coefficients in blocks of L, skipping Nyquist zeros at multiples of L.
        // Exploit symmetry to compute only the left half of the impulse response.
        // Assumes IMPULSE_RESPONSE_CENTER is divisible by INTERPOLATION_FACTOR.
        for (std::size_t block = 0; block < IMPULSE_RESPONSE_CENTER; block += L) {

            for (std::size_t k = block + 1; k < block + L; ++k) {

                const float n = static_cast<float>(k) - center;
                const float pi_n = M_PI*n;
                const float pi_n_inv = 1.0f/pi_n;
                const float ideal = std::sin(pi_n*R_inv)*pi_n_inv;
                const float norm  = n*center_inv;
                const float beta_sqrt_r = beta*std::sqrt(1.0f - norm*norm);
                const float win   = utils::float_arithmetic::besselI0(beta_sqrt_r)*I0_beta_inv;
                const float coeff = R*ideal*win;

                h[k] = coeff;
                h[FILTER_PROTOTYPE_LENGTH - k] = coeff; // symmetric
            }
        }

        // lim[n -> 0] sin(pi*n/R)/(pi*n) = 1/R.
        // The Kaiser window is exactly 1 at its center.
        h[IMPULSE_RESPONSE_CENTER] = R_inv;

        return h;
    }



    static float compute_kaiser_beta(const float atten) {

        if (atten >= 50.0) {
            return 0.1102*(atten - KAISER_SHAPE_PARAMETER_MATLAB_OFFSET);
        }
        if (atten > 21.0) {
            const float a = atten - 21.0f;
            return 0.5842*std::pow(a, 0.4) + 0.07886*a;
        }
        return 0.0;
    }




} // namespace channels::filter
