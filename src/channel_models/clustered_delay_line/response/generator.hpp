#pragma once

#include <array>

// Ray offset in 7.5 Table 7.5-3
constexpr std::array<float, 20> RAY_OFFSET_ALPHA = {
    0.0447f,  -0.0447f,   0.1413f,  -0.1413f,   0.2492f,  
   -0.2492f,   0.3715f,  -0.3715f,   0.5129f,  -0.5129f,
    0.6797f,  -0.6797f,   0.8844f,  -0.8844f,  -1.1481f,
   -1.1481f,   1.5195f,  -1.5195f,   2.1551f,  -2.1551f
};

namespace channels::cdl::response {

    //void generate_static_path_gains();
    //void generate_time_varying_path_gains();
    //void generate_ofdm_response();
    //void generate();
}