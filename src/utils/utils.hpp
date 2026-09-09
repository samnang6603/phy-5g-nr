#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <numbers>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include "../nr5g_common.hpp"

/***************** Constants ************************/
constexpr float KDEG2RAD = std::numbers::pi_v<float>/180.0f;
constexpr float DEG2RAD(float degrees) {
    return degrees*KDEG2RAD;
}
constexpr float LOG10_TOLERANCE = 300.0f;
constexpr float POW2DB(float x) {
    return (10.0f*std::log10(x) + LOG10_TOLERANCE) - LOG10_TOLERANCE;
}
constexpr float DB2POW(float x) {
    return std::pow(10.0f,x*0.1f);
}

namespace utils {

    namespace binary_arithmetic {

        std::vector<uint64_t> mapBin2Gray(const std::size_t M);

        uint16_t packBits(const uint8_t* bits, std::size_t bitsPerSymbol);

    }

    namespace files_import {

        template<typename T>
        std::vector<T> read_from_bin_files(const std::string& filename) {

            std::ifstream file(filename, std::ios::binary);

            if (!file) {
                throw std::runtime_error(
                    "Could not open file: " + filename);
            }

            // Read vector length
            uint64_t length = 0;

            file.read(reinterpret_cast<char*>(&length), sizeof(length));

            if (!file) {
                throw std::runtime_error(
                    "Failed to read vector length from: " + filename);
            }

            // Allocate vector
            std::vector<T> data(length);

            // Read vector data
            if (length > 0) {
                file.read(
                    reinterpret_cast<char*>(data.data()), 
                    static_cast<std::streamsize>(length*sizeof(T))
                );

                if (!file)
                {
                    throw std::runtime_error(
                        "Failed to read vector data from: " + filename);
                }
            }

            return data;
        }
    }

    namespace float_arithmetic {

        float calculateMeans(std::vector<float> x);

    }
}