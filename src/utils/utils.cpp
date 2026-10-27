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

        float calculateMeans(const std::vector<float> x) {

            float sum = 0.0f;
            for (const auto& v : x) {
                sum += v;
            }
            return sum/x.size();
        }

        float besselI0(const float x) {

            // Computes modified Bessel I_0 function of first kind
            // up to 20 finite terms (sufficient for beta around 6.75)

            // I_0(x) = sum_{k=0}^infinity ((x^2/4)^k)/(k!)^2
            // Let y = x^2/4 and T_k = y^k/(k!)^2
            // Since k! = k*(k-1)!, simplify to T_k = T_{k-1} * y/k^2
            // This recurrence avoids expensive factorial and pow() computations
            const std::size_t num_terms = BESSELI0_FLOAT64_TERM;
            const float y = 0.25f*x*x;
            float term = 1.0f; // T_0 = y^0/(0!)^2 = 1
            float sum = term;

            for (std::size_t i = 1; i < num_terms; ++i) {

                const float if32 = static_cast<float>(i);
                term *= y/(if32*if32); // T_i = T_{i-1} * y/i^2
                const float next = sum + term;
                sum = next;
            }

            return sum;
        }

        float dot(const float* a, const float* b, const std::size_t len) {

            float c = 0.0f;
            for (std::size_t i = 0; i < len; ++i) {
                c += (*a++)*(*b++);
            }
            return c;

        }

    }

}