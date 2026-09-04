#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>
#include "cblas.h"
#include "../../../nr5g_common.hpp"


static constexpr uint8_t BG1_LAYER_ROWS = 46;
static constexpr uint8_t BG2_LAYER_ROWS = 42;
static constexpr uint8_t BG1_MSG_CHUNKS = 22;
static constexpr uint8_t BG1_ALL_CHUNKS = 68;
static constexpr uint8_t BG1_PAR_CHUNKS = BG1_ALL_CHUNKS - BG1_MSG_CHUNKS;
static constexpr uint8_t BG2_MSG_CHUNKS = 10;
static constexpr uint8_t BG2_ALL_CHUNKS = 52;
static constexpr uint8_t BG2_PAR_CHUNKS = BG2_ALL_CHUNKS - BG2_MSG_CHUNKS;

namespace common::fec::ldpc {

    enum class BGN {BG1=1, BG2};

    class nrLDPC {
    public:
        nrLDPC(const uint64_t M, const float rate);
        ~nrLDPC() = default;

    private:
        uint64_t A_ = 0;       // information bit length
        BGN BGn_ = static_cast<BGN>(NULL);     // base graph number
        uint16_t Zc_ = 0;      // shift size
        uint8_t  set_idx_ = 0; // set index [0, 7]
        uint64_t K_ = 0;       // systematic bits length in base graph
        uint64_t N_ = 0;       // coded output length
        uint64_t n_F_ = 0;     // number of filler bits

        std::vector<uint16_t> edges_;  // shift connections from check nodes to variable nodes  
        std::vector<uint16_t> shifts_; // shifts for each correponding edge
        std::vector<uint16_t> layers_; // each layer is defined by edges which  the corresponding check node and  varible nodes) 
                                       // specified by start edge index and number of edges

    public:
        std::vector<uint8_t> encode(const std::vector<uint8_t>& msg) {

            // LDPC Encodes msg using ldpc_conf configurations

            uint8_t n_eqn, n_par;
            uint16_t n_msg_chunks, n_par_chunks;

            switch (BGn_) {
                case BGN::BG1:
                    n_eqn = BG1_LAYER_ROWS;
                    n_par = BG1_PAR_CHUNKS;
                    n_msg_chunks = static_cast<uint16_t>(BG1_MSG_CHUNKS);
                    n_par_chunks = static_cast<uint16_t>(BG1_PAR_CHUNKS);
                    break;
                
                case BGN::BG2:
                    n_eqn = BG2_LAYER_ROWS;
                    n_par = BG2_PAR_CHUNKS;
                    n_msg_chunks = static_cast<uint16_t>(BG2_MSG_CHUNKS);
                    n_par_chunks = static_cast<uint16_t>(BG2_PAR_CHUNKS);
                    break;

                default:
                    throw std::runtime_error("nrLDPC:encode >>> Invalid base graph");
            }

            // Validations
            const bool valid_msg_size = msg.size() == static_cast<std::size_t>(n_msg_chunks*Zc_);
            assert(valid_msg_size);
            const bool valid_eqn_par  = n_eqn == n_par_chunks;
            assert(valid_eqn_par);

            

            








            std::vector<uint8_t> coded_bits;



            return coded_bits;
        };
    };

    std::vector<uint8_t> ldpc_encode(
        const common::fec::ldpc::nrLDPC ldpc_conf, 
        const std::vector<uint8_t>& msg
    );
    
}