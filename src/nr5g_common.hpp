#pragma once

#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <complex>
#include <array>
#include <algorithm>
#include <fftw3.h>
#include <zlog.h>
#include <cassert>

namespace common {
    namespace modulation {
        enum class modulationType {BPSK, QPSK, PSK8, QAM, QAM16, QAM64, QAM128};
    }

    namespace fec {
        namespace ldpc {
            namespace encoder {
                
            
            }
            namespace decoder {
                enum class algorithm {BELIEF_PROPAGATION, MIN_SUM};
            }
        };
        namespace polar {

        }
    }
}

