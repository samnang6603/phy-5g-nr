#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include "../nr5g_common.hpp"

namespace utils {

    std::vector<uint64_t> mapBin2Gray(const std::size_t M);
    uint16_t packBits(const uint8_t* bits, std::size_t bitsPerSymbol);

}