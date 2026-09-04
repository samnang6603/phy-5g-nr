#pragma once

#include <vector>
#include "../../../nr5g_common.hpp"

/***************** Constants ************************/


/***************** Type Definitions ************/
namespace channels::cdl {

    struct Size {
        uint16_t M = 2;  // number of rows in antenna array
        uint16_t N = 2;  // number of cols in antenna array
        uint16_t P = 2;  // number of polarizations (1 or 2)
        uint16_t Mg = 1; // number of row of array panels
        uint16_t Ng = 1; // number of col of array panels
    };

    struct Spacing {
        float lambda_v = 0.5f; // vertical element spacing
        float lambda_h = 0.5f; // horizontal element spacing
        float dg_v = 1.0f; // vertical panel spacing
        float dg_h = 1.0f; // horizontal panel spacing
    };

    struct PolarizationAngles {
        float theta =  45.0f;
        float rho   = -45.0f;
    };

    struct ArrayOrientation {
        float alpha = 0.0f;
        float beta  = 0.0f;
        float gamma = 0.0f;
    };

    enum class ElementPattern { TR_38_901, ISOTROPIC };
    enum class PolarizationModel { MODEL1 = 1, MODEL2 };

    struct AntennaArray {
        Size struct_size;
        Spacing element_spacing;
        PolarizationAngles pol_angles;
        ArrayOrientation orientation;
        ElementPattern element = ElementPattern::TR_38_901;
        PolarizationModel pol_model = PolarizationModel::MODEL2;
    };

}