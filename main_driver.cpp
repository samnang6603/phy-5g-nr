#include <arm/types.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <zlog.h>
#include <cblas.h>
#include "simulations/channel_tests/cdl_tests/cdl_simulations.hpp"
#include "src/common/symbols_modulation/symbols_modulation.hpp"

#include "simulations/submodule_tests/ldpc_tests/ldpc_tests.hpp"


int main(void) {
    std::cout << "Hello World! " << std::tan(M_PI/4) << std::endl;

    std::vector<uint8_t> b{1,1,0,0,1,0,1,1};

    std::vector<std::complex<float>> bpsk_out(b.size(),0);

    std::size_t M = 4;

    std::vector<std::complex<float>> qpsk_out(b.size()/std::log2f(static_cast<float>(M)),0);

    common::modulation::modulateQPSK(qpsk_out, b, 0);

    //common::fec::ldpc::runSimulation1();

    channels::cdl::runSimulation1();


    std::cout << "Simulation Done!" << std::endl;


    return 0;
}