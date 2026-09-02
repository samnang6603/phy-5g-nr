#pragma once

#include <cstdint>
#include <vector>
#include "cblas.h"

namespace common::fec::ldpc {

    class nrLDPC {
    public:
        nrLDPC(uint64_t M, uint64_t rate);
        ~nrLDPC();

    private:
        uint64_t A_ = 0;       // information bit length
        uint64_t M_ = 0;       // codeword length
        uint8_t  BGn_ = 0;     // base graph number
        uint16_t Zc_ = 0;      // shift size
        uint8_t  set_idx_ = 0; // set index [0, 7]
        uint64_t K_ = 0;       // systematic bits length in base graph
        uint64_t N_ = 0;       // coded output length
        uint64_t n_F_ = 0;     // number of filler bits

        std::vector<uint16_t> edges_;
        std::vector<uint16_t> shifts_;
        std::vector<uint16_t> layers_;
    };

    std::vector<uint8_t> ldpc_encode(
        const common::fec::ldpc::nrLDPC ldpc_conf, 
        std::vector<uint8_t>& msg
    );
    
}