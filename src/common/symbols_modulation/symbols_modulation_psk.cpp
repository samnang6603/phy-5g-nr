#include "symbols_modulation.hpp"

namespace common {
    namespace modulation {

        template <std::size_t M>
        static void rotateConstellation(std::array<std::complex<float>, M>& constellation,
                                        const std::array<std::complex<float>, M>& baseConstellation,
                                        const float phaseOffset)
        {
            // Rotate constellations by phaseOffset amount

            const float cosPhase = std::cos(phaseOffset);
            const float sinPhase = std::sin(phaseOffset);

            for (std::size_t i = 0; i < M; ++i) {
                const float real = baseConstellation[i].real();
                const float imag = baseConstellation[i].imag();

                constellation[i].real(real*cosPhase - imag*sinPhase);
                constellation[i].imag(real*sinPhase + imag*cosPhase);
            }
        }

        void modulateBPSK(std::vector<std::complex<float>>& y,
                          const std::vector<uint8_t>& x)
        {
            // Modulate bits using BPSK

            y.resize(x.size());

            for (std::size_t i = 0; i < x.size(); ++i) {
                y[i].real(x[i] ? 1.0f : -1.0f);
                y[i].imag(0.0f);
            }
        }

        void modulateQPSK(std::vector<std::complex<float>>& y,
                          const std::vector<uint8_t>& x,
                          const float phaseOffset)
        {
            // Modulate symbols using QPSK

            constexpr std::size_t M = 4;
            constexpr std::size_t bitsPerSymbol = 2;

            if (x.size() % bitsPerSymbol != 0) {
                throw std::runtime_error("common::modulation::modulateQPSK(): Invalid number of input bits");
            }

            std::array<std::complex<float>, M> constellation;

            rotateConstellation(constellation, QPSK_CONSTELLATION_LUT, phaseOffset);

            y.resize(x.size()/bitsPerSymbol);

            for (std::size_t i = 0; i < y.size(); ++i) {
                const uint16_t sym = utils::binary_arithmetic::packBits(&x[i*bitsPerSymbol], bitsPerSymbol);
                y[i] = constellation[sym];
            }
        }

        void modulate8PSK(std::vector<std::complex<float>>& y,
                          const std::vector<uint8_t>& x,
                          const float phaseOffset)
        {
            // Modulate symbols using 8-PSK

            constexpr std::size_t M = 8;
            constexpr std::size_t bitsPerSymbol = 3;

            if (x.size() % bitsPerSymbol != 0) {
                throw std::runtime_error("common::modulation::modulate8PSK(): Invalid number of input bits");
            }

            std::array<std::complex<float>, M> constellation;

            rotateConstellation(constellation, PSK8_CONSTELLATION_LUT, phaseOffset);

            y.resize(x.size()/bitsPerSymbol);

            for (std::size_t i = 0; i < y.size(); ++i) {
                const uint16_t sym = utils::binary_arithmetic::packBits(&x[i*bitsPerSymbol], bitsPerSymbol);
                y[i] = constellation[sym];
            }
        }
    }
}