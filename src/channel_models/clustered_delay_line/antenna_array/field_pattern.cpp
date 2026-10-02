#include <array>
#include <cmath>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cblas.h>
#include "../../../utils/utils.hpp"
#include "antenna.hpp"

namespace channels::cdl::antenna::field_pattern {

/***************** Constants ************************/
static constexpr float TR38_901_SLA_V = 30.0f;          // Side-Lobe Attenuation (dB)
static constexpr float TR38_901_THETA_3dB = 65.0f;      // Vertical Half-Power Beam-Width (degrees)
static constexpr float TR38_901_THETA_3dB_INV = 1.0f/TR38_901_THETA_3dB;
static constexpr float TR38_901_A_M = 30.0f;            // Front-to-back attenuation ratio (dB)
static constexpr float TR38_901_PHI_3dB = 65.0f;        // Horizontal Half-Power Beam-Width (degrees)
static constexpr float TR38_901_PHI_3dB_INV = 1.0f/TR38_901_PHI_3dB;
static constexpr float TR38_901_G_MAX = 8.0f;           // Maximum directional gain of an element (dBi)
static constexpr float PSI_DEGENERATE_THRESH = 1E-5f;   // Degenerate tolerance protection
static constexpr float TWO_PI = M_PI*2.0f;

}

/***************** Function Implementations ************/
namespace channels::cdl::antenna::field_pattern {

    template<ElementPatternOption PowMode, PolarizationModelOption PolMode>
    static std::vector<float> compute_pattern_kernel(
        const std::vector<float>& theta_p,
        const std::vector<float>& phi_p,
        const float cos_zeta,
        const float sin_zeta
    );

    static std::vector<float> compute_polarization_field_pattern(
        const ElementPatternOption pow_mode,
        const PolarizationModelOption pol_mode,
        const std::vector<float>& theta_p,
        const std::vector<float>& phi_p,
        const float zeta
    );

    namespace nlos {

        void compute_field_term(
            float* field_term,
            const antenna::AntennaArrayConfig& ant_array_conf,
            const float theta,
            const float phi,
            const std::size_t ant_idx
        ) {
            // Compute field term of NLoS clusters in GCS
            // Using Equation 7.1-11. The process of this function involves finding psi
            // and F_prime from polarization and field pattern

            /*
            NOTE in TR 38.901:  
            - theta and phi are Azimuth and Zenith, respectively, in GCS
            - theta_prime and phi_prime as Azimuth and Zenith, respectively, in LCS
            This distinction is critical for transforming field vectors and applying 
            element pattern gains.

            Orientation angles of the array
            Bearing and Downtilt angles are the combination of array and element
            orientation
            Bearing angle (x-axis)
            */
            const geometry::ArrayOrientationConfig& init_ort_conf = ant_array_conf.InitialOrientation;
            const geometry::ArrayOrientationConfig& ort_conf = ant_array_conf.Orientation;
            
            // Compute rotation matrix from LCS to GCS for the *initial* orientation.
            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> rot_init = 
                geometry::compute_lcs2gcs(init_ort_conf);

            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> rot = 
                geometry::compute_lcs2gcs(ort_conf);

            /*
            Transform the target orientation into the local frame of the initial 
            orientation.
            That is: R = RInit' * ROriented * RInit
            This performs a change of basis:
            - Maps ROriented (in GCS) into the LCS of the reference orientation.
            - The result, R, expresses the target orientation *as seen from* the 
                reference LCS.
            
            Extra: To transform a matrix A into another basis B, we use
            A_in_B = B^T * A * B
            
            WHY go back LCS frame?
            BECAUSE: We care about how the rotated antenna (element + array) looks 
            from the antenna’s starting orientation — because field vectors, 
            element patterns, and mutual coupling are defined in the antenna's own 
            local frame, not the world frame.
            */
            const float* A = rot.data();
            const float* B = rot_init.data();
            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> R = 
                geometry::multiply_3x3_matrices_BtAB(A, B);

            // sph_unit: [rho, theta, phi] <=> [radius, azimuth, zenith]

        }
    } // namespace nlos

    namespace los {
        
        void compute_field_term(
            float* field_term,
            const antenna::AntennaArrayConfig& ant_array_conf,
            const geometry::los::SphericalDirection sph_dir_conf,
            const std::size_t ant_idx
        ) {

            // Compute field term of LoS clusters in GCS
            // Using Equation 7.1-11. The process of this function involves finding psi
            // and F_prime from polarization and field pattern

            /*
            NOTE in TR 38.901:  
            - theta and phi are Azimuth and Zenith, respectively, in GCS
            - theta_prime and phi_prime as Azimuth and Zenith, respectively, in LCS
            This distinction is critical for transforming field vectors and applying 
            element pattern gains.

            Orientation angles of the array
            Bearing and Downtilt angles are the combination of array and element
            orientation
            Bearing angle (x-axis)
            */
            const geometry::ArrayOrientationConfig& init_ort_conf = ant_array_conf.InitialOrientation;
            const geometry::ArrayOrientationConfig& ort_conf = ant_array_conf.Orientation;
            
            // Compute rotation matrix from LCS to GCS for the *initial* orientation.
            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> rot_init = 
                geometry::compute_lcs2gcs(init_ort_conf);

            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> rot = 
                geometry::compute_lcs2gcs(ort_conf);

            /*
            Transform the target orientation into the local frame of the initial 
            orientation.
            That is: R = RInit' * ROriented * RInit
            This performs a change of basis:
            - Maps ROriented (in GCS) into the LCS of the reference orientation.
            - The result, R, expresses the target orientation *as seen from* the 
                reference LCS.
            
            Extra: To transform a matrix A into another basis B, we use
            A_in_B = B^T * A * B
            
            WHY go back LCS frame?
            BECAUSE: We care about how the rotated antenna (element + array) looks 
            from the antenna’s starting orientation — because field vectors, 
            element patterns, and mutual coupling are defined in the antenna's own 
            local frame, not the world frame.
            */
            const float* A = rot.data();
            const float* B = rot_init.data();
            const std::array<float, geometry::NUM_ELEMENT_3D_GRID> R = 
                geometry::multiply_3x3_matrices_BtAB(A, B);

            // sph_unit: [rho, theta, phi] <=> [radius, azimuth, zenith]
            //const std::array<float, geometry::NUM_3D_AXIS> rho_hat = 
            //    geometry::get_LoS_spherical_unit_vector(phi, theta);

            // Compute theta_prime and phi_prime according to TR 38.901 Eqn 7.1-7 and 7.1-8
            const std::array<float, geometry::NUM_3D_AXIS>& rhat = sph_dir_conf.rhat;
            float theta_prime = std::acos(rhat[2]);
            float phi_prime   = std::atan2(rhat[1], rhat[0]);
            if (theta_prime == 0.0f) {
                phi_prime = 0.0f;
            }

            // Compute psi, the angular displacement between two pairs of
            // unit vectors, according to Equation 7.1-12
            // First, we need to find the unit vector of theta, phi, and theta_p.
            // These are the Cartesian representation of the spherical unit vectors
            const float costheta = sph_dir_conf.costheta;
            const float sintheta = sph_dir_conf.sintheta;

            const float cosphi = sph_dir_conf.cosphi;
            const float sinphi = sph_dir_conf.sinphi;

            const float cosphi_prime = std::cos(phi_prime);
            const float sinphi_prime = std::sin(phi_prime);

            const float costheta_prime = std::cos(theta_prime);
            const float sintheta_prime = std::sin(theta_prime);

            // The unit vectors theta_hat, phi_hat, and theta_prime_hat
            const std::array<float, geometry::NUM_3D_AXIS> theta_hat = {
                costheta*cosphi,
                costheta*sinphi, 
                -sintheta
            };

            const std::array<float, geometry::NUM_3D_AXIS> phi_hat = {
                -sinphi, 
                cosphi, 
                0.0f
            };

            const std::array<float, geometry::NUM_3D_AXIS> theta_prime_hat = {
                costheta_prime*cosphi_prime,
                costheta_prime*sinphi_prime,
                -sintheta_prime
            };

            // Using Equation 7.1-12, simplified by hand to get the least computational
            // expression
            const float tmp0 = R[0]*theta_prime_hat[0] + R[3]*theta_prime_hat[1] + R[6]*theta_prime_hat[2];
            const float tmp1 = R[1]*theta_prime_hat[0] + R[4]*theta_prime_hat[1] + R[7]*theta_prime_hat[2];
            const float tmp2 = R[2]*theta_prime_hat[0] + R[5]*theta_prime_hat[1] + R[8]*theta_prime_hat[2];

            const float y = phi_hat[0]*tmp0 + phi_hat[1]*tmp1 + phi_hat[2]*tmp2;
            const float x = theta_hat[0]*tmp0 + theta_hat[1]*tmp1 + theta_hat[2]*tmp2;

            const float psi = std::atan2(y, x);

            const PatternConfig& pattern_conf = ant_array_conf.FieldPattern.FieldEffect;
            const PolarizationAnglesConfig& pol_conf = ant_array_conf.PolarizationAngles;

            const std::vector<float> theta_prime_vec = {RAD2DEG(theta_prime)};
            const std::vector<float> phi_prime_vec   = {RAD2DEG(phi_prime)};

            // Element's antenna index polarization angle selector
            const float zeta = pol_conf.PolarizationOrientationMap[geometry::NUM_3D_AXIS*ant_idx + 2];

            // Compute field pattern
            std::vector<float> F = compute_polarization_field_pattern(
                pattern_conf.Element, 
                pattern_conf.PolarizationModel,
                theta_prime_vec,
                phi_prime_vec,
                zeta
            );

            // Finally, we have everything we need to compute fieldTerm in GCS,
            // using equation 7.1-11
            const float cospsi = std::cos(psi);
            const float sinpsi = std::sin(psi);

            //              [cospsi  -sinpsi]
            // fieldTerm1 = |               | * F
            //              [sinpsi   cospsi]
            // Alternatively, simplied form, hand derived for optimization

            field_term[0] = F[0]*cospsi - F[1]*sinpsi;
            field_term[1] = F[0]*sinpsi + F[1]*cospsi;
        }

        void get_location_term(
            std::complex<float>* loc_term,
            const std::array<float, geometry::NUM_3D_AXIS>& rhat,
            const std::vector<float>& dbar,
            const float lambda_0_inv,
            const std::size_t ant_idx
        ) {
            // Get location term

            const float v1 = rhat[0]*dbar[ant_idx*geometry::NUM_3D_AXIS];
            const float v2 = rhat[1]*dbar[ant_idx*geometry::NUM_3D_AXIS + 1];
            const float v3 = rhat[2]*dbar[ant_idx*geometry::NUM_3D_AXIS + 2];
            const float sumv = TWO_PI*(v1 + v2 + v3)*lambda_0_inv;
            const std::complex<float> c(std::cos(sumv), std::sin(sumv));

            std::fill(loc_term, loc_term + NUMBER_OF_RAYS_PER_CLUSTER, c);
        }
    } // namespace los


    float get_polarization_angle(
        std::size_t antennaIdx,
        const AntennaArrayConfig& ant_array_conf
    ) {
        // Returns the polarization slant angle for a given antenna-port index 
        // based on the configured polarization ordering

        const geometry::SizeConfig& s = ant_array_conf.Size;
        const field_pattern::PolarizationAnglesConfig& pol = ant_array_conf.PolarizationAngles;

        const std::size_t elements_per_pol =
            static_cast<std::size_t>(s.M)*s.N;

        const std::size_t p = (antennaIdx/elements_per_pol) % s.P;

        return (p == 0) ? pol.theta : pol.rho;
    }

    static std::vector<float> compute_polarization_field_pattern(
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

        std::vector<float> F; // Field pattern output

        const float zeta_r = DEG2RAD(zeta);
        const float cos_zeta = std::cos(zeta_r);
        const float sin_zeta = std::sin(zeta_r);

        switch (pol_mode) {
            case PolarizationModelOption::MODEL1:
                switch (pow_mode) {
                    case ElementPatternOption::TR_38_901:
                        F = compute_pattern_kernel<
                            ElementPatternOption::TR_38_901,
                            PolarizationModelOption::MODEL1
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta
                        );
                        break;

                    case ElementPatternOption::ISOTROPIC:
                        F = compute_pattern_kernel<
                            ElementPatternOption::ISOTROPIC,
                            PolarizationModelOption::MODEL1
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta
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
                        F = compute_pattern_kernel<
                            ElementPatternOption::TR_38_901,
                            PolarizationModelOption::MODEL2
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta
                        );
                        break;

                    case ElementPatternOption::ISOTROPIC:
                        F = compute_pattern_kernel<
                            ElementPatternOption::ISOTROPIC,
                            PolarizationModelOption::MODEL2
                        >(
                            theta_p,
                            phi_p,
                            cos_zeta,
                            sin_zeta
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

    template<ElementPatternOption PowMode>
    static float compute_power_pattern(
        const float theta_p_deg,
        const float phi_p_deg
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
            const float tmp0 = (theta_p_deg - 90)*TR38_901_THETA_3dB_INV;
            const float tmp1 = 12.0f*tmp0*tmp0;
            const float A_EV = -std::min(tmp1, TR38_901_SLA_V);

            // Antenna element horizontal radiation pattern (dB)
            const float tmp2 = phi_p_deg*TR38_901_PHI_3dB_INV;
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

    template<ElementPatternOption PowMode, PolarizationModelOption PolMode>
    static std::vector<float> compute_pattern_kernel(
        const std::vector<float>& theta_p,
        const std::vector<float>& phi_p,
        const float cos_zeta,
        const float sin_zeta
    ) {
        // Compute field pattern kernel

        static_assert(
            PolMode == PolarizationModelOption::MODEL1 ||
            PolMode == PolarizationModelOption::MODEL2,
            "compute_polarization >>> Invalid polarization model"
        );

        const std::size_t theta_len = theta_p.size();

        std::vector<float> F(NUM_MAX_POLARIZATION);

        for (std::size_t i = 0; i < theta_len; ++i) {

            const float prad = compute_power_pattern<PowMode>(theta_p[i], phi_p[i]);
            const float prad_sqrt = std::sqrt(prad);

            if constexpr (PolMode == PolarizationModelOption::MODEL1) {
                // TR 38.901 7.3.2 Model-1

                // Rotation matrix elements cos(Psi) and sin(Psi) for an angular
                // displacement of Psi due to the orientation of the LCS w.r.t. the GCS
                // See Equation 7.3-3 cos(phi) and sin(phi)
                const float theta_p_r = DEG2RAD(theta_p[i]);
                const float sin_theta = std::sin(theta_p_r);
                const float cos_theta = std::cos(theta_p_r);

                const float phi_p_r = DEG2RAD(phi_p[i]);
                const float sin_phi = std::sin(phi_p_r);
                const float cos_phi = std::cos(phi_p_r);

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
                F[i] = prad_sqrt*cos_psi;
                F[i + theta_len] = prad_sqrt*sin_psi;
            }

            if constexpr (PolMode == PolarizationModelOption::MODEL2) {

                // TR 38.901 Equation 7.3-4
                F[i] = prad_sqrt*cos_zeta;

                // TR 38.901 Equation 7.3-5
                F[i + theta_len] = prad_sqrt*sin_zeta;
            }
        }

        return F;
    }
    

} // namespace channels::cdl::antenna::field_pattern

