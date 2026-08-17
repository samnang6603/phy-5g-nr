#pragma once

#include <algorithm>
#include <vector>

namespace common {
    
    namespace fec {

        namespace ldpc {

            namespace gf2 {
                
                struct GF2MatrixConfig {
                    uint32_t rows = 0;              // destLen
                    uint32_t cols = 0;              // length(source)
                    std::vector<uint32_t> colPtr;   // size cols + 1
                    std::vector<uint16_t> rowInd;
                    std::size_t destWords() const { 
                        return (rows + 63) / 64; 
                    }  // for bit-packed dest
                };

                template<typename T>
                GF2MatrixConfig createMatrixGF2(
                    const T* rowInd,
                    const T* rowloc,
                    const T* colsum,
                    const uint32_t cols,
                    const uint32_t rows,
                    bool dedupCols = false
                );

            } // namespace gf2
        } // namespace ldpc
    } // namespace fec
} // namespace common

