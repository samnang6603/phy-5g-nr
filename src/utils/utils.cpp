#include "utils.hpp"
#include <stdexcept>

namespace utils {

    namespace binary_arithmetic {

        std::vector<uint64_t> mapBin2Gray(const std::size_t M) {

            std::vector<uint64_t> index(M,0);

            for (std::size_t i = 0; i < M; ++i) {
                index[i ^ (i >> 1)] = i;
            }

            return index;
        }

        uint16_t packBits(const uint8_t* bits, std::size_t bitsPerSymbol)
        {
            // pack input bits into symbols

            uint16_t sym = 0;

            for (std::size_t i = 0; i < bitsPerSymbol; ++i) {
                sym = static_cast<uint16_t>((sym << 1) | bits[i]);
            }

            return sym;
        }
    }

    namespace float_arithmetic {

        float calculateMeans(std::vector<float> x) {

            float sum;
            for (const auto& v : x) {
                sum += v;
            }

            return sum/x.size();

        }

    }

}