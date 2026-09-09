#include <stdexcept>
#include <vector>
#include <algorithm>
#include "antenna_structure.hpp"
#include "../cdl.hpp"

/***************** Constants ************************/
static constexpr float TR38_901_SLA_V = 30.0f;     /// Side-Lobe Attenuation (dB)
static constexpr float TR38_901_THETA_3dB = 65.0f; // Vertical Half-Power Beam-Width (degrees)
static constexpr float TR38_901_THETA_3dB_INV = 1/TR38_901_THETA_3dB;
static constexpr float TR38_901_A_M = 30.0f;       // Front-to-back attenuation ratio (dB)
static constexpr float TR38_901_PHI_3dB = 65.0f;   // Horizontal Half-Power Beam-Width (degrees)
static constexpr float TR38_901_PHI_3dB_INV = 1/TR38_901_PHI_3dB; 
static constexpr float TR38_901_G_MAX = 8.0f;      // Maximum directional gain of an element (dBi)

/***************** Function Implementations ************/
namespace channels::cdl::antenna::field_effect {

    void compute_power_pattern(
        const ElementPattern powmode,
        const std::vector<float> theta_p,
        const std::vector<float> phi_p
    ) {
        // Compute A_prime, antenna gain pattern Antenna element radiation 
        // pattern is described in TR 38.901 Section 7.3 table 7.3-1

        const std::size_t theta_len = theta_p.size();

        std::vector<float> A(theta_len);
        std::vector<float> A_EV(theta_len);
        std::vector<float> A_EH(theta_len);
        
        switch (powmode) {
            case ElementPattern::TR_38_901:
            case ElementPattern::BS:
                for (std::size_t i = 0; i < theta_len; ++i) {

                    // Antenna element vertical radiation pattern (dB)
                    const float tmp = theta_p[i]*TR38_901_THETA_3dB_INV;
                    const float tmp1 = 12*tmp*tmp;
                    A_EV[i] = -std::min(tmp1, TR38_901_SLA_V);

                    // Antenna element horizontal radiation pattern (dB)
                    const float tmp2 = phi_p[i]*TR38_901_THETA_3dB_INV;
                    const float tmp3 = 12*tmp2*tmp2;
                    A_EH[i] = -std::min(tmp3,TR38_901_A_M);

                    // Combining method for 3D antenna element pattern (dB)
                    const float tmp4 = -std::min(-(A_EV[i] + A_EH[i]), TR38_901_A_M);

                    // Incorporate maximum gain and convert to linear power
                    A[i] = std::pow(10,(tmp4 + TR38_901_G_MAX)*0.1f);
                }
                break;
            
            case ElementPattern::ISOTROPIC:
                // Equally radiate power in all directions, pattern is circle
                std::fill(A.begin(),A.end(),1.0f);
                break;

            default:
                throw std::runtime_error("compute_power_pattern >>> Invalid power pattern");
        }

    }

}