#include <string>
#include "ldpc_tests.hpp"
#include "../../../src/utils/utils.hpp"
#include "../../../src/common/fec/ldpc/ldpc.hpp"


namespace common::fec::ldpc {

    void runSimulation1(void) {

        const std::string filename = "input_bin_files/ldpc_msg_seed777.bin";
        std::vector<uint8_t> msg = utils::files_import::read_from_bin_files<uint8_t>(filename);

        const uint64_t M = 4000;
        const float rate = 0.5;

        nrLDPC ldpc_conf(M,rate);
        
    }
   
}