#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>

/***************** Constants ************************/
constexpr int FILTER_DELAY = 7;
// 16 Provides best balance between accuracy, 
// passband characteristics, and computational efficiency
constexpr uint16_t FILTER_LENGTH = 16;
constexpr uint16_t FILTER_LENGTH_HALF = FILTER_LENGTH/2;
constexpr uint16_t INTERPOLATION_FACTOR = 50;
constexpr uint16_t FILTER_PROTOTYPE_LENGTH = INTERPOLATION_FACTOR*FILTER_LENGTH;
constexpr uint16_t IMPULSE_RESPONSE_CENTER = FILTER_PROTOTYPE_LENGTH/2;
constexpr float STOPBAND_ATTENUATION = 70.0f; // dB
constexpr float KAISER_SHAPE_PARAMETER_MATLAB_OFFSET = 8.71f;

// Fractional-delay candidates ordered to match MATLAB's polyphase filter layout.
// Phase 0 represents zero delay; subsequent phases descend from (L-1)/L to 1/L
// because the sinc peak shifts left as the polyphase index increases.
// The final 1.0 handles rounding near the next integer-sample boundary,
// ensuring a maximum fractional-delay quantization error of 1/(2*L) = 0.01.
constexpr std::array<float, INTERPOLATION_FACTOR + 1> INTERPOLATION_RANGE = {
    0.0f   , 0.9800f, 0.9600f, 0.9400f, 0.9200f, 0.9000f, 
    0.8800f, 0.8600f, 0.8400f, 0.8200f, 0.8000f, 0.7800f, 
    0.7600f, 0.7400f, 0.7200f, 0.7000f, 0.6800f, 0.6600f, 
    0.6400f, 0.6200f, 0.6000f, 0.5800f, 0.5600f, 0.5400f, 
    0.5200f, 0.5000f, 0.4800f, 0.4600f, 0.4400f, 0.4200f, 
    0.4000f, 0.3800f, 0.3600f, 0.3400f, 0.3200f, 0.3000f, 
    0.2800f, 0.2600f, 0.2400f, 0.2200f, 0.2000f, 0.1800f, 
    0.1600f, 0.1400f, 0.1200f, 0.1000f, 0.0800f, 0.0600f, 
    0.0400f, 0.0200f, 1.0000f
};


namespace channels::filter {

    struct ChannelFilterConfig {
        float fs = 15360000.0f;
        std::array<float, FILTER_PROTOTYPE_LENGTH> channel_filter_coeffs;
        std::vector<float> path_filter_coeffs;
        std::array<float, FILTER_PROTOTYPE_LENGTH> initial_conditions;
        std::array<float, FILTER_PROTOTYPE_LENGTH> final_conditions;
    };

    void design_channel_filter(
        ChannelFilterConfig& chan_filt_conf,
        const std::vector<float> path_delays
    );




    


} // namespace channels::filter