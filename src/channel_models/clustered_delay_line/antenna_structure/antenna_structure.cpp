
#include <vector>
#include <cmath>
#include <array>
#include <cstddef>
#include "cblas.h"
#include "../../../utils/utils.hpp"
#include "antenna_structure.hpp"

/***************** Constants ************************/


/***************** Implementations ************************/
namespace channels::cdl::antenna {

    static void arrange_layout(AntennaArrayConfig& ant_array);
    static std::array<float, NUM_ELEMENT_3D_GRID> compute_lcs2gcs_zbroadside(const ArrayOrientationConfig& ort);
    static std::vector<float> compute_reoriented_antenna_positions(
        const std::array<float, NUM_ELEMENT_3D_GRID>& V,
        const std::vector<float>& pos,
        const SizeConfig& s
    );

    AntennaSystemConfig initializeAntennaStructure(
        AntennaArrayConfig& tx_ant_array,
        AntennaArrayConfig& rx_ant_array
    ) {

        // Arrange antenna structure and create spatial antenna position map 

        arrange_layout(tx_ant_array);
        arrange_layout(rx_ant_array);

        AntennaSystemConfig ant_sys_conf;
        ant_sys_conf.TransmitAntennaArray = tx_ant_array;
        ant_sys_conf.ReceiveAntennaArray  = rx_ant_array;

        return ant_sys_conf;

    }

    float get_polarization_angle(
        std::size_t antennaIdx,
        const AntennaArrayConfig& ant_array
    ) {
        // Returns the polarization slant angle for a given antenna-port index 
        // based on the configured polarization ordering

        const SizeConfig& s = ant_array.ArraySize;
        const PolarizationAnglesConfig& pol = ant_array.PolarizationAngles;

        const std::size_t elements_per_pol =
            static_cast<std::size_t>(s.M)*s.N;

        const std::size_t p = (antennaIdx/elements_per_pol) % s.P;

        return (p == 0) ? pol.theta : pol.rho;
    }

    static void arrange_layout(AntennaArrayConfig& ant_array) {

        // Arrange antenna spatial layout

        const auto& s = ant_array.ArraySize;
        const auto& spc = ant_array.ElementSpacing;
        const auto& ort = ant_array.Orientation;
        const auto& pol = ant_array.PolarizationAngles;

        std::vector<float> spacing = {
            spc.d_v,
            spc.d_h,
            spc.dg_v,
            spc.dg_h,
        };

        std::vector<uint16_t> size = {
            s.M,
            s.N,
            s.Mg,
            s.Ng,
            s.P,
        };

        std::vector<float> pos(s.num_antenna_ports()*NUM_3D_AXIS,0.0f);

        // The pos array is arranged in this way:
        // [ [x0, y0, z0, x1, y1, z1, ..., xM, yM, zM],
        //   [x0, y0, z0, x1, y1, z1, ..., xN, yN, zN],
        //   [x0, y0, z0, x1, y1, z1, ..., xMg, yMg, zMg],
        //   [x0, y0, z0, x1, y1, z1, ..., xNg, yNg, zNg] ]
        // positional axis offset refers to the coordinate position
        // starting with x-axis for M (row), pos_axis_offset = 0
        // for y-axis for N (col), pos_axis_offset = 1
        // Note that z-axis is not accounted for
        std::size_t pos_axis_offset = 0; 

        std::size_t spc_idx = 0;
        std::size_t pos_idx = 0;
        for (std::size_t i = 0; i < size.size() - 1; ++i) {

            const auto& v = size[i];

            float accum_tmp = spacing[spc_idx++];
            std::vector<float> tmp;
            tmp.push_back(accum_tmp);

            for (std::size_t b = 0; b < v - 1; ++b) {
                accum_tmp += accum_tmp;
                tmp.push_back(accum_tmp);
            }

            float mean_tmp = utils::float_arithmetic::calculateMeans(tmp);

            for (std::size_t b = 0; b < v; ++b) {
                pos[pos_idx] = tmp[b] - mean_tmp;
                pos_idx += NUM_3D_AXIS;
            }

            // After M (pos_axis_offset = 0), then N (pos_axis_offset = 1), reset to
            // pos_axis_offset = 0 for Mg and then pos_axis_offset = 1 for Ng
            pos_axis_offset = (pos_axis_offset + 1) % 2;
            pos_idx += pos_axis_offset;

        }

        std::array<float, NUM_ELEMENT_3D_GRID> vgcs_broadside = compute_lcs2gcs_zbroadside(ort);
        ant_array.ElementPositions = compute_reoriented_antenna_positions(vgcs_broadside, pos, s);
        
        // TODO: polarization ormap needed?

    }

    static std::array<float, NUM_ELEMENT_3D_GRID> compute_lcs2gcs_zbroadside(const ArrayOrientationConfig& ort) {

        const float a = ort.alpha;
        const float b = ort.beta;
        const float g = ort.gamma;

        const float ca = std::cos(DEG2RAD(a));
        const float cb = std::cos(DEG2RAD(b));
        const float cg = std::cos(DEG2RAD(g));

        const float sa = std::sin(DEG2RAD(a));
        const float sb = std::sin(DEG2RAD(b));
        const float sg = std::sin(DEG2RAD(g));

        std::array<float, NUM_ELEMENT_3D_GRID> r;

        // Column 0 = -(original column 2)
        r[0] = -(ca*sb*cg + sa*sg);
        r[1] = -(sa*sb*cg - ca*sg);
        r[2] = -(cb*cg);

        // Column 1 = original column 1
        r[3] =  ca*sb*sg - sa*cg;
        r[4] =  sa*sb*sg + ca*cg;
        r[5] =  cb*sg;

        // Column 2 = original column 0
        r[6] =  ca*cb;
        r[7] =  sa*cb;
        r[8] = -sb;

        return r;
    }

    static inline std::vector<float> transform_pos(
        const std::vector<float>& pos,
        const std::array<float, NUM_ELEMENT_3D_GRID>& V
    ) {
        // Applies a 3×3 coordinate transformation matrix to each 3D position
        // vector in the flat position array.

        std::vector<float> out(pos.size());

        for (std::size_t i = 0; i < pos.size(); i += NUM_3D_AXIS) {
            const float r0 = pos[i + 0];
            const float r1 = pos[i + 1];
            const float r2 = pos[i + 2];

            out[i + 0] = V[0]*r0 + V[3]*r1 + V[6]*r2;
            out[i + 1] = V[1]*r0 + V[4]*r1 + V[7]*r2;
            out[i + 2] = V[2]*r0 + V[5]*r1 + V[8]*r2;
        }

        return out;
    }

    static std::vector<float> compute_reoriented_antenna_positions(
        const std::array<float, NUM_ELEMENT_3D_GRID>& V,
        const std::vector<float>& pos,
        const SizeConfig& s
    ) {

        // Computes the 3D position of every antenna port by combining the per-dimension position 
        // offsets after applying the global coordinate transformation.

        std::vector<float> pos_r = transform_pos(pos, V);

        const std::size_t M  = s.M;
        const std::size_t N  = s.N;
        const std::size_t P  = s.P;
        const std::size_t Mg = s.Mg;
        const std::size_t Ng = s.Ng;

        const std::size_t offM  = 0;
        const std::size_t offN  = offM  + NUM_3D_AXIS*M;
        const std::size_t offP  = offN  + NUM_3D_AXIS*N;
        const std::size_t offMg = offP  + NUM_3D_AXIS*P;
        const std::size_t offNg = offMg + NUM_3D_AXIS*Mg;

        std::vector<float> reoriented(NUM_3D_AXIS*s.num_antenna_ports());

        std::size_t k = 0;

        for (std::size_t ng = 0; ng < Ng; ++ng) {
            const std::size_t iNg = offNg + NUM_3D_AXIS*ng;

            for (std::size_t mg = 0; mg < Mg; ++mg) {
                const std::size_t iMg = offMg + NUM_3D_AXIS*mg;

                const float mgng0 = pos_r[iMg + 0] + pos_r[iNg + 0];
                const float mgng1 = pos_r[iMg + 1] + pos_r[iNg + 1];
                const float mgng2 = pos_r[iMg + 2] + pos_r[iNg + 2];

                for (std::size_t p = 0; p < P; ++p) {
                    const std::size_t iP = offP + NUM_3D_AXIS*p;

                    const float pmgng0 = pos_r[iP + 0] + mgng0;
                    const float pmgng1 = pos_r[iP + 1] + mgng1;
                    const float pmgng2 = pos_r[iP + 2] + mgng2;

                    for (std::size_t n = 0; n < N; ++n) {
                        const std::size_t iN = offN + NUM_3D_AXIS*n;

                        const float npmgng0 = pos_r[iN + 0] + pmgng0;
                        const float npmgng1 = pos_r[iN + 1] + pmgng1;
                        const float npmgng2 = pos_r[iN + 2] + pmgng2;

                        for (std::size_t m = 0; m < M; ++m) {
                            const std::size_t iM = offM + NUM_3D_AXIS*m;
                            const std::size_t k3 = NUM_3D_AXIS*k;

                            reoriented[k3 + 0] = pos_r[iM + 0] + npmgng0;
                            reoriented[k3 + 1] = pos_r[iM + 1] + npmgng1;
                            reoriented[k3 + 2] = pos_r[iM + 2] + npmgng2;

                            ++k;
                        }
                    }
                }
            }
        }

        return reoriented;
    }

    



}