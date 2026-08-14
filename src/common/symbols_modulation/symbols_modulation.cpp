#include "symbols_modulation.hpp"
#include <stdexcept>

namespace common {
    namespace modulation {

        /************************* Private static functions ***********************/

        void modulate(
            std::vector<std::complex<float>>& y,
            const std::vector<uint8_t>& x,
            const common::modulation::modulationType modType,
            const float phaseOffset
        ) {
            switch (modType) {
                case common::modulation::modulationType::BPSK:
                    modulateBPSK(y,x);
                    break;

                case common::modulation::modulationType::QPSK:
                    modulateQPSK(y,x,phaseOffset);
                    break;

                case common::modulation::modulationType::PSK8:
                    modulate8PSK(y,x,phaseOffset);
                    break;

                /*
                case common::modulation::modulationType::QAM16:
                    modulate16QAM(...);
                    break;
                */
                default:
                    throw std::runtime_error("common::modulation::modulate(): Invalid modulation");
                
        
            }
        }
    }
}