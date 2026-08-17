#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include <span>
#include "psk_lut.hpp"
#include "../../nr5g_common.hpp"
#include "../../utils/utils.hpp"


namespace common {
    namespace modulation {

        void modulate(
            std::vector<std::complex<float>>& y,
            const std::vector<uint8_t>& x,
            const common::modulation::modulationType modType
        );

        void modulateBPSK(std::vector<std::complex<float>>& y, const std::vector<uint8_t>& x);
        void modulateQPSK(std::vector<std::complex<float>>& y, const std::vector<uint8_t>& x, const float phaseOffset);
        void modulate8PSK(std::vector<std::complex<float>>& y, const std::vector<uint8_t>& x, const float phaseOffset);        


    } // namespace modulation
} // namespace common