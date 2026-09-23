#pragma once

#include <array>
#include <cstdint>

constexpr std::uint32_t MT19937_STATE_SIZE = 624U;
constexpr std::uint32_t MT19937_PERIOD_OFFSET = 397U;
constexpr std::uint32_t MT19937_TWIST_OFFSET = MT19937_STATE_SIZE - MT19937_PERIOD_OFFSET;

constexpr std::uint32_t MT19937_MATRIX_A = 0x9908B0DFU;
constexpr std::uint32_t MT19937_UPPER_MASK = 0x80000000U;
constexpr std::uint32_t MT19937_LOWER_MASK = 0x7FFFFFFFU;

constexpr std::uint32_t MT19937_SEED_MULTIPLIER = 1812433253U;

constexpr std::uint32_t MT19937_TEMPERING_MASK_B = 0x9D2C5680U;
constexpr std::uint32_t MT19937_TEMPERING_MASK_C = 0xEFC60000U;

constexpr std::uint32_t MATLAB_DEFAULT_MT19937_SEED = 5489U;

constexpr double MT19937_DOUBLE_MULTIPLIER = 67108864.0;       // 2^26
constexpr double MT19937_DOUBLE_DIVISOR = 9007199254740992.0; // 2^53


class MATLABMT19937
{
public:
    explicit MATLABMT19937(std::uint32_t seed) {
        // MATLAB RandStream('mt19937ar','Seed',0)
        // uses the MT19937 default seed 5489.
        const std::uint32_t actualSeed =
            (seed == 0U) ? MATLAB_DEFAULT_MT19937_SEED : seed;

        state_[0] = actualSeed;

        for (std::uint32_t i = 1; i < MT19937_STATE_SIZE; ++i) {
            state_[i] =
                MT19937_SEED_MULTIPLIER *
                (state_[i - 1] ^ (state_[i - 1] >> 30)) +
                i;
        }

        index_ = MT19937_STATE_SIZE;
    }

    std::uint32_t nextUInt32() {
        if (index_ >= MT19937_STATE_SIZE) {
            twist();
        }

        std::uint32_t y = state_[index_++];

        y ^= (y >> 11);
        y ^= (y << 7) & MT19937_TEMPERING_MASK_B;
        y ^= (y << 15) & MT19937_TEMPERING_MASK_C;
        y ^= (y >> 18);

        return y;
    }

    double rand() {

        double value;

        do {
            // MATLAB mt19937ar full-precision rand uses
            // two 32-bit MT19937 outputs to construct
            // one 53-bit double.

            const std::uint32_t a = nextUInt32() >> 5;
            const std::uint32_t b = nextUInt32() >> 6;

            value = (
                    static_cast<double>(a) *
                    MT19937_DOUBLE_MULTIPLIER +
                    static_cast<double>(b)
                )/MT19937_DOUBLE_DIVISOR;

            // MATLAB rand returns values strictly inside (0,1).
            // If MT19937 generates exactly zero, MATLAB discards
            // it and generates another value.
        } while (value == 0.0);
        
        return value;
    }

private:
    void twist() {
        std::uint32_t y;

        for (std::uint32_t i = 0; i < MT19937_TWIST_OFFSET; ++i)
        {
            y =
                (state_[i] & MT19937_UPPER_MASK) |
                (state_[i + 1U] & MT19937_LOWER_MASK);

            state_[i] =
                state_[i + MT19937_PERIOD_OFFSET] ^
                (y >> 1U) ^
                ((y & 1U) ? MT19937_MATRIX_A : 0U);
        }

        for (
            std::uint32_t i = MT19937_TWIST_OFFSET;
            i < MT19937_STATE_SIZE - 1U;
            ++i
        ) {
            y =
                (state_[i] & MT19937_UPPER_MASK) |
                (state_[i + 1U] & MT19937_LOWER_MASK);

            state_[i] =
                state_[i - MT19937_TWIST_OFFSET] ^
                (y >> 1U) ^
                ((y & 1U) ? MT19937_MATRIX_A : 0U);
        }

        y =
            (state_[MT19937_STATE_SIZE - 1U] & MT19937_UPPER_MASK) |
            (state_[0] & MT19937_LOWER_MASK);

        state_[MT19937_STATE_SIZE - 1U] =
            state_[MT19937_PERIOD_OFFSET - 1U] ^
            (y >> 1U) ^
            ((y & 1U) ? MT19937_MATRIX_A : 0U);

        index_ = 0U;
    }

    std::array<std::uint32_t, MT19937_STATE_SIZE> state_{};
    std::uint32_t index_ = MT19937_STATE_SIZE;
};