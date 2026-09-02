#include "ldpc.hpp"

namespace common::fec::ldpc {

    nrLDPC::nrLDPC(uint64_t M, uint64_t rate) {
        
    }

    std::vector<uint8_t> ldpc_encode(
        const common::fec::ldpc::nrLDPC ldpc_conf, 
        const std::vector<uint8_t>& msg
    ) {
        // LDPC Encodes msg using ldpc_conf configurations



    }

}