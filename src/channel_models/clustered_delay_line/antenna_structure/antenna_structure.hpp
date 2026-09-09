#pragma once

#include <array>
#include <vector>
#include "../../../nr5g_common.hpp"
#include "../../../utils/utils.hpp"

/***************** Constants ************************/
static constexpr uint8_t NUM_3D_AXIS = 3;
static constexpr uint8_t NUM_ELEMENT_3D_GRID = 9;
static constexpr uint8_t PHYSICAL_ANTENNA_ARRAY_PROPERTIES = 5;

/***************** Type Definitions ************/
namespace channels::cdl::antenna {

    struct SizeConfig {
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

    struct SpacingConfig {
        float d_v = 0.5f; // vertical element spacing
        float d_h = 0.5f; // horizontal element spacing
        float dg_v = 1.0f; // vertical panel spacing
        float dg_h = 1.0f; // horizontal panel spacing
    };

    struct PolarizationAnglesConfig {
        float theta =  45.0f; // degrees
        float rho   = -45.0f; // degrees

        float operator[](std::size_t p) const noexcept {
            return (p == 0) ? theta : rho;
        }
    };

    struct ArrayOrientationConfig {
        // All in degrees
        float alpha = 0.0f; 
        float beta  = 0.0f;
        float gamma = 0.0f;
    };

    enum class ElementPatternOption { 
        TR_38_901,
        ISOTROPIC 
    };
    
    enum class PolarizationModelOption { 
        MODEL1 = 1, 
        MODEL2
    };

    struct AntennaFieldEffectConfig {
        ElementPatternOption Element = ElementPatternOption::TR_38_901;
        PolarizationModelOption PolarizationModel = PolarizationModelOption::MODEL2;
        std::vector<float> F;
    };

    struct AntennaArrayConfig {
        SizeConfig ArraySize;
        SpacingConfig ElementSpacing;
        PolarizationAnglesConfig PolarizationAngles;
        ArrayOrientationConfig Orientation; 
        AntennaFieldEffectConfig FieldEffect;
        std::vector<float> ElementPositions;
    };

    struct AntennaSystemConfig {
        std::size_t NumInputSignals = 8;
        std::size_t NumOutputSignals = 2;
        AntennaArrayConfig TransmitAntennaArray;
        AntennaArrayConfig ReceiveAntennaArray;
    };
        
}

/***************** Function Declarations ************/
namespace channels::cdl::antenna {

    AntennaSystemConfig initializeAntennaStructure(
        AntennaArrayConfig& tx_ant_array,
        AntennaArrayConfig& rx_ant_array
    );

    float get_polarization_angle(
        std::size_t antennaIdx,
        const AntennaArrayConfig& ant_array
    );

    namespace field_effect {
        
        std::vector<float> compute_polarization_field_pattern(
            const ElementPatternOption pow_mode,
            const PolarizationModelOption pol_mode,
            const std::vector<float>& theta_p,
            const std::vector<float>& phi_p,
            const float zeta
        );
    }
}