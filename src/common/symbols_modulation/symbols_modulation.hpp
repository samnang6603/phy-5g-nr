#include <complex>
#include <vector>

namespace common {
    namespace modulation {

        void bpsk(std::vector<std::complex<float>>& y, const std::vector<bool>& x);
        void qpsk(std::vector<std::complex<float>>& y, const std::vector<bool>& x);
        void qam(std::vector<std::complex<float>>& y, const std::vector<bool>& x);
    }
}