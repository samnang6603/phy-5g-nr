#include "gf2.hpp"
#include <cstdint>

namespace common {
    
    namespace fec {

        namespace ldpc {

            namespace gf2 {

                unsigned ctz64(const uint64_t x) {

                    // Count trailing zeros
                    return static_cast<unsigned>(__builtin_ctzll(x));
                }

                template<typename T>
                GF2MatrixConfig createMatrixGF2(
                    const T* rowInd,
                    const T* rowloc,
                    const T* colsum,
                    const uint32_t cols,
                    const uint32_t rows,
                    bool dedupCols = false
                ) {

                    // Matrix multiplication over finite Galois Field 2 GF(2)

                    GF2MatrixConfig gf2mat_conf;
                    gf2mat_conf.cols = cols;
                    gf2mat_conf.rows = rows;
                    gf2mat_conf.colPtr.resize(cols + 1);
                    std::size_t nnz = 0;

                    for (uint32_t c = 0; c < cols; ++c) {
                        // Calculate total number of non-zeros
                        nnz += reinterpret_cast<std::size_t>(colsum);
                    }

                    // Reserve gf2 matrix row indices
                    gf2mat_conf.rowInd.reserve(nnz);

                    std::vector<uint16_t> tmp;
                    for (uint32_t c = 0; c < cols; ++c) {
                        gf2mat_conf.colPtr[c] = static_cast<uint32_t>(gf2mat_conf.rowInd.size());
                        const T* p = rowInd + static_cast<size_t>(rowloc[c]);
                        const T  n = colsum[c];

                        if (!dedupCols) {
                            for (T b = 0; b < n; ++b) {
                                gf2mat_conf.rowInd.push_back(static_cast<uint16_t>(p[b]));
                            }
                        } else {
                            tmp.assign(p, p + n);
                            std::sort(tmp.begin(), tmp.end());
                            tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
                            gf2mat_conf.rowInd.insert(gf2mat_conf.rowInd.end(), tmp.begin(), tmp.end());
                        }
                    }
                    gf2mat_conf.colPtr[cols] = static_cast<uint32_t>(gf2mat_conf.rowInd.size());
                    return gf2mat_conf;
                }

                void multiplyMatrixNaive(
                    uint8_t* __restrict dest,
                    const GF2MatrixConfig& matrix_conf,
                    const uint8_t* __restrict src
                ) {
                    // Computes the modulo-2 matrix product of a vector (specified by source) and a matrix
                    
                    std::memset(dest, 0, matrix_conf.rows);
                    for (uint32_t c = 0; c < matrix_conf.cols; ++c) {
                        if (!src[c]) {
                            continue;
                        }
                        const uint32_t start = matrix_conf.colPtr[c];
                        const uint32_t end = matrix_conf.colPtr[c + 1];
                        for (uint32_t k = start; k < end; ++k) {
                            dest[matrix_conf.rowInd[k]] ^= 1U;
                        }
                    }
                }

                void multiplyMatrixBitPacked(
                    uint64_t* __restrict dest,
                    const GF2MatrixConfig& matrix_conf,
                    const uint8_t* __restrict src
                ) {
                    // Computes the modulo-2 matrix product of a vector (specified by source) and a matrix
                    
                    std::memset(dest, 0, matrix_conf.destWords() * sizeof(uint64_t));
                    for (uint32_t c = 0; c < matrix_conf.cols; ++c) {
                        if (!src[c]) { 
                            continue;
                        }
                        const uint32_t beg = matrix_conf.colPtr[c], end = matrix_conf.colPtr[c + 1];
                        for (uint32_t k = beg; k < end; ++k) {
                            const uint32_t r = matrix_conf.rowInd[k];
                            dest[r >> 6] ^= 1ULL << (r & 63);
                        }
                    }
                }

                void substituteEquationMatrix(
                    uint8_t* __restrict x,
                    const GF2MatrixConfig& matrix_conf,
                    const int8_t direction
                ) {
                    // Solves a system of linear equations by forward/backward substitution (if source == dest).

                    if (direction > 0) {
                        for (uint32_t c = 0; c < matrix_conf.cols; ++c) {
                            if (!x[c]) continue;
                            for (uint32_t k = matrix_conf.colPtr[c]; k < matrix_conf.colPtr[c+1]; ++k)
                                x[matrix_conf.rowInd[k]] ^= 1u;
                        }
                    } else {
                        for (uint32_t c = matrix_conf.cols; c-- > 0; ) {
                            if (!x[c]) continue;
                            for (uint32_t k = matrix_conf.colPtr[c]; k < matrix_conf.colPtr[c+1]; ++k)
                                x[matrix_conf.rowInd[k]] ^= 1u;
                        }
                    }
                }
            } // namespace gf2
        } // namespace ldpc
    } // namespace fec
} // namespace common