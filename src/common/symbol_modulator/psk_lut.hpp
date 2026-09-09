#pragma once

#include <array>
#include <complex>

constexpr float INV_SQRT2 = 0.7071067811865475244f;

/*
 * LUTs are indexed by the packed binary value of the input bits.
 *
 * QPSK Gray mapping around the unit circle:
 *
 *   00 ->   0 deg
 *   01 ->  90 deg
 *   11 -> 180 deg
 *   10 -> 270 deg
 *
 * Since packBits() returns binary indices:
 *
 *   00 -> 0
 *   01 -> 1
 *   10 -> 2
 *   11 -> 3
 *
 * the LUT order is:
 *
 *   index 0 -> 00 ->   0 deg
 *   index 1 -> 01 ->  90 deg
 *   index 2 -> 10 -> 270 deg
 *   index 3 -> 11 -> 180 deg
 */
constexpr std::array<std::complex<float>, 4> QPSK_CONSTELLATION_LUT = {
    {
        { 1.0f,  0.0f},
        { 0.0f,  1.0f},
        { 0.0f, -1.0f},
        {-1.0f,  0.0f}
    }
};

/*
 * 8-PSK Gray mapping around the unit circle:
 *
 *   000 ->   0 deg
 *   001 ->  45 deg
 *   011 ->  90 deg
 *   010 -> 135 deg
 *   110 -> 180 deg
 *   111 -> 225 deg
 *   101 -> 270 deg
 *   100 -> 315 deg
 *
 * LUT is indexed by packed binary input:
 *
 *   000 -> 0
 *   001 -> 1
 *   010 -> 2
 *   011 -> 3
 *   100 -> 4
 *   101 -> 5
 *   110 -> 6
 *   111 -> 7
 */
constexpr std::array<std::complex<float>, 8> PSK8_CONSTELLATION_LUT = {
    {
        { 1.0f,      0.0f     },   // 000
        { INV_SQRT2,  INV_SQRT2 },   // 001
        {-INV_SQRT2,  INV_SQRT2 },   // 010
        { 0.0f,      1.0f     },   // 011
        { INV_SQRT2, -INV_SQRT2 },   // 100
        { 0.0f,     -1.0f     },   // 101
        {-1.0f,      0.0f     },   // 110
        {-INV_SQRT2, -INV_SQRT2 }    // 111
    }
};