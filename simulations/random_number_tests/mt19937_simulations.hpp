#pragma once

#include <array>
#include <iostream>
#include <iomanip>
#include "../../src/channel_models/matlab_mt19937.hpp"

inline void runRandomSimulation1(std::uint32_t seed) {

    MATLABMT19937 randomStream(seed);

    std::array<double, 10> values{};

    for (double& value : values)
    {
        value = randomStream.rand();
    }

    std::cout << std::setprecision(17);

    for (const double value : values)
    {
        std::cout << value << '\n';
    }

}