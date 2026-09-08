
#include <vector>
#include <cmath>
#include "../../../utils/utils.hpp"
#include "antenna_structure.hpp"

/***************** Constants ************************/
static constexpr float DEG2RAD(float degrees) {
    return degrees*M_PI/180.0f;
}

/***************** Implementations ************************/
namespace channels::cdl::antenna {

    static void arrange_layout(const AntennaArray& ant_array);
    static std::vector<float> compute_lcs2gcs_zbroadside(const ArrayOrientation& ort);

    void arrangeStructure(
        AntennaArray& ant_array
    ) {

        arrange_layout(ant_array);

    }

    static void arrange_layout(const AntennaArray& ant_array) {

        // Arrange antenna layout

        const auto& siz = ant_array.struct_size;
        const auto& spc = ant_array.element_spacing;
        const auto& ort = ant_array.orientation;

        std::vector<float> spacing = {
            spc.d_v,
            spc.d_h,
            spc.dg_v,
            spc.dg_h,
        };

        std::vector<uint16_t> size = {
            siz.M,
            siz.N,
            siz.Mg,
            siz.Ng,
            siz.P,
        };

        std::vector<float> pos(siz.num_antenna_ports*NUM_3D_AXIS,0.0f);

        // The pos array is arranged in this way:
        // [ [x0, y0, z0, x1, y1, z1, ..., xM, yM, zM],
        //   [x0, y0, z0, x1, y1, z1, ..., xN, yN, zN],
        //   [x0, y0, z0, x1, y1, z1, ..., xMg, yMg, zMg],
        //   [x0, y0, z0, x1, y1, z1, ..., xNg, yNg, zNg]
        // ]
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

        static std::vector<float> vgcs_broadside = compute_lcs2gcs_zbroadside(ort);

        

        



    }

    static std::vector<float> compute_lcs2gcs_zbroadside(const ArrayOrientation& ort) {

        const float a = ort.alpha;
        const float b = ort.beta;
        const float g = ort.gamma;

        const float ca = std::cos(DEG2RAD(a));
        const float cb = std::cos(DEG2RAD(b));
        const float cg = std::cos(DEG2RAD(g));

        const float sa = std::sin(DEG2RAD(a));
        const float sb = std::sin(DEG2RAD(b));
        const float sg = std::sin(DEG2RAD(g));

        std::vector<float> r(NUM_ELEMENT_3D_GRID);

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

    



}