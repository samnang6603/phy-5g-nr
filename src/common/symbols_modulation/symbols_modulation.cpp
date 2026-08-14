#include "symbols_modulation.hpp"

namespace common {
    namespace modulation {

        void bpsk(std::vector<std::complex<float>>& y, const std::vector<bool>& x) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                y[i].real(x[i] ? 1.0f : -1.0f);
            }
        }

        void qpsk(std::vector<std::complex<float>>& y, const std::vector<bool>& x) {

        }

        void qam(std::vector<std::complex<float>>& y, const std::vector<bool>& x) {

        }
    }
}