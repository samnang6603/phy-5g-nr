#pragma once

#include <cstdint>
#include <vector>
#include "cblas.h"

namespace common::fec::ldpc {

    enum BGN {BG1=1, BG2};

    class nrLDPC {
    public:
        nrLDPC(const uint64_t M, const float rate);
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

        std::vector<uint16_t> edges_;  // shift connections from check nodes to variable nodes  
        std::vector<uint16_t> shifts_; // shifts for each correponding edge
        std::vector<uint16_t> layers_; // each layer is defined by edges which  the corresponding check node and  varible nodes) 
                                       // specified by start edge index and number of edges
    };

    std::vector<uint8_t> ldpc_encode(
        const common::fec::ldpc::nrLDPC ldpc_conf, 
        const std::vector<uint8_t>& msg
    );
    
}