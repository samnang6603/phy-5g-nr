#include <vector>
#include <array>
#include "../../nr5g_common.hpp"

#pragma once

/***************** Constants ************************/
constexpr uint8_t SYMBOLS_PER_SLOT = 14;

namespace modules {

    class Carrier {
    public:
        std::size_t NumSubcarriers = 624;
        std::size_t NumFFT = 1024;
        float SampleRate = 15360000.0f;
        float SubcarrierSpacing = 15; // KHz
        std::size_t NumSizeGrid = 52;
        std::size_t NSlot = 0;
        std::size_t NFrame = 0;

        uint16_t SymbolsPerSlot() const noexcept {
            return symbols_per_slot_;
        }

        uint16_t SlotsPerSubframe() const noexcept {
            return slots_per_subframe_;
        }

        uint16_t SlotsPerFrame() const noexcept {
            return slots_per_frame_;
        }

    private:
        uint16_t symbols_per_slot_ = 14;
        uint16_t slots_per_subframe_ = 1;
        uint16_t slots_per_frame_ = 10;
    };

} // namespace modules
