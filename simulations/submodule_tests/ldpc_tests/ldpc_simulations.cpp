#include "ldpc_tests.hpp"
#include "../../../src/utils/utils.hpp"
#include <string>

namespace common::fec::ldpc {

    void runSimulation1() {

        const std::string filename = "input_bin_files/ldpc_msg_seed777.bin";
        std::vector<uint8_t> msg = utils::files_import::read_from_bin_files<uint8_t>(filename);

        

    }
   
}