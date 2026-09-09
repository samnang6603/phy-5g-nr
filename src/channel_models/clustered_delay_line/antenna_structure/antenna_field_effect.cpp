#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include "../../../utils/utils.hpp"
#include "../cdl.hpp"
#include "antenna_structure.hpp"

/***************** Constants ************************/
static constexpr float TR38_901_SLA_V = 30.0f;          // Side-Lobe Attenuation (dB)
static constexpr float TR38_901_THETA_3dB = 65.0f;      // Vertical Half-Power Beam-Width (degrees)
static constexpr float TR38_901_THETA_3dB_INV = 1.0f/TR38_901_THETA_3dB;
static constexpr float TR38_901_A_M = 30.0f;            // Front-to-back attenuation ratio (dB)
static constexpr float TR38_901_PHI_3dB = 65.0f;        // Horizontal Half-Power Beam-Width (degrees)
static constexpr float TR38_901_PHI_3dB_INV = 1.0f/TR38_901_PHI_3dB;
static constexpr float TR38_901_G_MAX = 8.0f;           // Maximum directional gain of an element (dBi)
static constexpr float PSI_DEGENERATE_THRESH = 1E-5f;   // Degenerate tolerance protection

/***************** Function Implementations ************/
namespace channels::cdl::antenna::field_effect {

    template<ElementPatternOption PowMode>
    static float compute_power_pattern(
        const float theta_p,
        const float phi_p
    ) {
        // Compute A_prime, antenna gain pattern Antenna element radiation
        // pattern is described in TR 38.901 Section 7.3 table 7.3-1

        static_assert(
            PowMode == ElementPatternOption::TR_38_901 ||
            PowMode == ElementPatternOption::ISOTROPIC,
            "compute_power_pattern >>> Invalid power pattern"
        );

        if constexpr (
            PowMode == ElementPatternOption::TR_38_901
        ) {
            // Antenna element vertical radiation pattern (dB)
            const float tmp0 = theta_p*TR38_901_THETA_3dB_INV;
            const float tmp1 = 12.0f*tmp0*tmp0;
            const float A_EV = -std::min(tmp1, TR38_901_SLA_V);

            // Antenna element horizontal radiation pattern (dB)
            const float tmp2 = phi_p*TR38_901_PHI_3dB_INV;
            const float tmp3 = 12.0f*tmp2*tmp2;
            const float A_EH = -std::min(tmp3, TR38_901_A_M);

            // Combining method for 3D antenna element pattern (dB)
            const float tmp4 = -std::min(-(A_EV + A_EH), TR38_901_A_M);

            // Incorporate maximum gain and convert to linear power
            return std::pow(10.0f, (tmp4 + TR38_901_G_MAX)*0.1f);
        }

        if constexpr (PowMode == ElementPatternOption::ISOTROPIC) {
            // Equally radiate power in all directions, pattern is circle
            return 1.0f;
        }
    }

    template<
        ElementPatternOption PowMode,
        PolarizationModelOption PolMode
    >
    static void compute_pattern_kernel(
        const std::vector<float>& theta_p,
        const std::vector<float>& phi_p,
        const float cos_zeta,
        const float sin_zeta,
        std::vector<float>& F
    ) {
        static_assert(
            PolMode == PolarizationModelOption::MODEL1 ||
            PolMode == PolarizationModelOption::MODEL2,
            "compute_polarization >>> Invalid polarization model"
        );

        const std::size_t theta_len = theta_p.size();

        for (std::size_t i = 0; i < theta_len; ++i) {

            const float prad = compute_power_pattern<PowMode>(theta_p[i], phi_p[i]);

            if constexpr (PolMode == PolarizationModelOption::MODEL1) {
                // TR 38.901 7.3.2 Model-1

                // Rotation matrix elements cos(Psi) and sin(Psi) for an angular
                // displacement of Psi due to the orientation of the LCS w.r.t. the GCS
                // See Equation 7.3-3 cos(phi) and sin(phi)
                const float theta_p_d = DEG2RAD(theta_p[i]);
                const float sin_theta = std::sin(theta_p_d);
                const float cos_theta = std::cos(theta_p_d);

                const float phi_p_d = DEG2RAD(phi_p[i]);
                const float sin_phi = std::sin(phi_p_d);
                const float cos_phi = std::cos(phi_p_d);

                const float tmp = sin_zeta*sin_phi;
                const float tmp0 = cos_zeta*cos_theta;
                const float tmp1 = tmp*sin_theta;
                const float tmp2 = (tmp0 - tmp1);

                const float denom = std::sqrt(1.0f - tmp2*tmp2);
                const float denom_inv = 1.0f/denom;

                float cos_psi = (cos_zeta*sin_theta + tmp*cos_theta)*denom_inv;
                float sin_psi = (sin_zeta*cos_phi)*denom_inv;

                // Assume vertical polarization in cases where the transformation
                // degenerates, this is, when theta_p is near zeta or 180-zeta and
                // phi_p is near -90 or 90, respectively. When zeta = 0/180, any
                // value of phi_p makes the transformation degenerate. The threshold
                // (10^-5) is the upper bound of the magnitude error outside a region of
                // radius ~10^-4 degrees around the singularity. Within that region, the
                // polarization angle error can be arbitrary.
                const float hypot = std::sqrt(cos_psi*cos_psi + sin_psi*sin_psi);

                const bool degen =
                    std::isnan(cos_psi) ||
                    std::isnan(sin_psi) ||
                    (std::fabs(hypot - 1.0f) > PSI_DEGENERATE_THRESH);

                if (degen) {
                    cos_psi = 1.0f;
                    sin_psi = 0.0f;
                }

                // Equation 7.3-3 evaluated assuming F_phi_pp = 0
                /*
                    [F_theta_p]     [cos(phi) -sin(phi)]   [F_theta_pp]
                    |         |  =  |                  | * |          |
                    [ F_phi_p ]     [sin(phi)  cos(phi)]   [ F_phi_pp ]
                */
                F[i] = prad*cos_psi;
                F[i + theta_len] = prad*sin_psi;
            }

            if constexpr (PolMode == PolarizationModelOption::MODEL2) {

                // TR 38.901 Equation 7.3-4
                F[i] = prad*cos_zeta;

                // TR 38.901 Equation 7.3-5
                F[i + theta_len] = prad*sin_zeta;
            }
        }
    }

    std::vector<float> compute_polarization_field_pattern(
        const ElementPatternOption pow_mode,
        const PolarizationModelOption pol_mode,
        const std::vector<float>& theta_p,
        const std::vector<float>& phi_p,
        const float zeta
    ) {
        // Compute field pattern of elements in LCS

        if (theta_p.size() != phi_p.size()) {
            throw std::invalid_argument(
                "compute_polarization_field_pattern >>> "
                "theta_p and phi_p must have the same size"
            );
        }

        const std::size_t theta_len = theta_p.size();

        std::vector<float> F(theta_len*2); // Field pattern output

        const float zeta_d = DEG2RAD(zeta);
        const float cos_zeta = std::cos(zeta_d);
        const float sin_zeta = std::sin(zeta_d);

        switch (pol_mode) {
            case PolarizationModelOption::MODEL1:
                switch (pow_mode) {
                    case ElementPatternOption::TR_38_901:
                        compute_pattern_kernel<
                            ElementPatternOption::TR_38_901,
                            PolarizationModelOption::MODEL1
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta,
                            F
                        );
                        break;

                    case ElementPatternOption::ISOTROPIC:
                        compute_pattern_kernel<
                            ElementPatternOption::ISOTROPIC,
                            PolarizationModelOption::MODEL1
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta,
                            F
                        );
                        break;

                    default:
                        throw std::invalid_argument(
                            "compute_polarization_field_pattern >>> Invalid power pattern"
                        );
                }
                break;

            case PolarizationModelOption::MODEL2:
                switch (pow_mode) {
                    case ElementPatternOption::TR_38_901:
                        compute_pattern_kernel<
                            ElementPatternOption::TR_38_901,
                            PolarizationModelOption::MODEL2
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta,
                            F
                        );
                        break;

                    case ElementPatternOption::ISOTROPIC:
                        compute_pattern_kernel<
                            ElementPatternOption::ISOTROPIC,
                            PolarizationModelOption::MODEL2
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta,
                            F
                        );
                        break;

                    default:
                        throw std::invalid_argument(
                            "compute_power_pattern >>> Invalid power pattern"
                        );
                }
                break;

            default:
                throw std::invalid_argument(
                    "compute_polarization >>> Invalid polarization model"
                );
        }

        return F;
    }

}