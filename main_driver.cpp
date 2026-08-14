#include <arm/types.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <zlog.h>
#include <cblas.h>
#include "src/common/symbols_modulation/symbols_modulation.hpp"


int main(void) {
    std::cout << "Hello World! " << std::tan(M_PI/4) << std::endl;

    std::vector<bool> v;
    v.push_back(true);
    v.push_back(false);
    v.push_back(true);
    v.push_back(false);
    v.push_back((false));

    std::vector<std::complex<float>> bpsk_out(5,0);

    common::modulation::bpsk(bpsk_out, v);


    return 0;
}