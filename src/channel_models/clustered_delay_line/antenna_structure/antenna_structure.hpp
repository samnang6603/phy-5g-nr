#pragma once

#include <array>
#include <vector>
#include "../../../nr5g_common.hpp"

/***************** Constants ************************/
static constexpr uint8_t NUM_3D_AXIS = 3;
static constexpr uint8_t NUM_ELEMENT_3D_GRID = 9;
static constexpr uint8_t PHYSICAL_ANTENNA_ARRAY_PROPERTIES = 5;

/***************** Type Definitions ************/
namespace channels::cdl::antenna {

    struct Size {
        uint16_t M = 2;  // number of rows in antenna array
        uint16_t N = 2;  // number of cols in antenna array
        uint16_t P = 2;  // number of polarizations (1 or 2)
        uint16_t Mg = 1; // number of row of array panels
        uint16_t Ng = 1; // number of col of array panels

        std::size_t num_antennas() const noexcept {
            return static_cast<std::size_t>(M)*N*Mg*Ng;
        }

        std::size_t num_antenna_ports() const noexcept {
            return num_antennas()*P;
        }

        std::size_t num_spatial_values() const noexcept {
            return NUM_3D_AXIS*num_antenna_ports();
        }
    };

    struct Spacing {
        float d_v = 0.5f; // vertical element spacing
        float d_h = 0.5f; // horizontal element spacing
        float dg_v = 1.0f; // vertical panel spacing
        float dg_h = 1.0f; // horizontal panel spacing
    };

    struct PolarizationAngles {
        // All in degrees
        float theta =  45.0f;
        float rho   = -45.0f;

        float operator[](std::size_t p) const noexcept {
            return (p == 0) ? theta : rho;
        }
    };

    struct ArrayOrientation {
        // All in degrees
        float alpha = 0.0f; 
        float beta  = 0.0f;
        float gamma = 0.0f;
    };

    enum class ElementPattern { 
        TR_38_901,
        BS,
        ISOTROPIC 
    };
    
    enum class PolarizationModel { 
        MODEL1 = 1, 
        MODEL2
    };

    struct AntennaArray {
        Size struct_size;
        Spacing element_spacing;
        PolarizationAngles pol_angles;
        ArrayOrientation orientation;
        ElementPattern element = ElementPattern::TR_38_901;
        PolarizationModel pol_model = PolarizationModel::MODEL2;
    };

    
}

/***************** Function Declarations ************/
namespace channels::cdl::antenna {

    void arrangeStructure(
        AntennaArray& ant_array
    );

    float get_polarization_angle(
        std::size_t antennaIdx,
        const AntennaArray& ant_array
    );

    
}