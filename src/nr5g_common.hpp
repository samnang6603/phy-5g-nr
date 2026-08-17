#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <array>
#include <fftw3.h>
#include <zlog.h>

namespace common {
    namespace modulation {
        enum class modulationType {BPSK, QPSK, PSK8, QAM, QAM16, QAM64, QAM128};
    }

    namespace fec {
        namespace ldpc {
            namespace encoder {
                enum class baseGraph {BG1, BG2};
            }
            namespace decoder {
                enum class algorithm {BELIEF_PROPAGATION, MIN_SUM};
            }
        };
        namespace polar {

        }
    }
}

