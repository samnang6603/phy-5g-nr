#pragma once

#include <ranges>
#include <array>
#include <vector>
#include "../../channel_models_common.hpp"
#include "../../../utils/utils.hpp"


/***************** Type Definitions ************/
namespace channels::cdl::antenna {

    inline constexpr std::size_t NUM_MAX_POLARIZATION = 2;

    namespace geometry {

        /***************** Constants ************************/
        inline constexpr uint8_t NUM_3D_AXIS = 3;
        inline constexpr uint8_t NUM_ELEMENT_3D_GRID = 9;
        inline constexpr uint8_t PHYSICAL_ANTENNA_ARRAY_PROPERTIES = 5;

        // Fixed LoS transform of target orientation into the LCS frame of the initial orientation
        constexpr std::array<float, NUM_ELEMENT_3D_GRID> LOS_LCS_INIT_ROTATION = {
            1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f 
        }; 

        /***************** Type Definitions ************************/
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

        struct PositionConfig {
            // Initialize with (0, 0, 0)
            float x = 0.0f;
            float y = 0.0f;
            float z = 0.0f;
        };

        struct SpacingConfig {
            float d_v = 0.5f; // vertical element spacing
            float d_h = 0.5f; // horizontal element spacing
            float dg_v = 1.0f; // vertical panel spacing
            float dg_h = 1.0f; // horizontal panel spacing
        };

        struct ArrayOrientationConfig {
            // All in degrees
            float alpha = 0.0f; // bearing (x-axis)
            float beta  = 0.0f; // downtilt (z-axis)
            float gamma = 0.0f; // slant (y-axis)
        };

        namespace los {

            struct SphericalDirection {
                std::array<float, NUM_3D_AXIS> rhat;

                float sintheta;
                float costheta;
                float sinphi;
                float cosphi;
            };
        } // namespace los

        namespace nlos {

            struct SphericalDirections {
                std::vector<float> rhat;

                std::vector<float> sintheta;
                std::vector<float> costheta;
                std::vector<float> sinphi;
                std::vector<float> cosphi;
            };
        } // namespace nlos

    } // namespace geometry

    namespace field_pattern {

        struct PolarizationAnglesConfig {
            float theta =  45.0f; // degrees
            float rho   = -45.0f; // degrees

            float operator[](std::size_t p) const noexcept {
                return (p == 0) ? theta : rho;
            }

            std::vector<float> PolarizationOrientationMap;
        };

        enum class ElementPatternOption { 
            TR_38_901,
            ISOTROPIC 
        };
        
        enum class PolarizationModelOption { 
            MODEL1 = 1, 
            MODEL2
        };

        struct PatternConfig {
            ElementPatternOption Element = ElementPatternOption::TR_38_901;
            PolarizationModelOption PolarizationModel = PolarizationModelOption::MODEL2;
            std::vector<float> F;
        };

        struct Config {
            PatternConfig FieldEffect;
            std::vector<float> ElementPositions;
        };

    } // namespace field_pattern

    struct AntennaArrayConfig {
        geometry::SizeConfig Size;
        geometry::PositionConfig Position;
        geometry::SpacingConfig ElementSpacing;
        geometry::ArrayOrientationConfig InitialOrientation{0.0f,0.0f,0.0f};
        geometry::ArrayOrientationConfig Orientation;
        field_pattern::PolarizationAnglesConfig PolarizationAngles;
        field_pattern::Config FieldPattern;
    };

    struct AntennaSystemConfig {
        std::size_t NumInputSignals;
        std::size_t NumOutputSignals;
        AntennaArrayConfig TransmitAntennaArray;
        AntennaArrayConfig ReceiveAntennaArray;
        float CarrierFrequency = 4E+9f;
        
        float Wavelength() const noexcept {
            return SPEED_OF_LIGHT / CarrierFrequency;
        }
    };
        
} // namespace antenna

/***************** Function Declarations ************/
namespace channels::cdl::antenna {

    namespace geometry {

        void initialize(AntennaSystemConfig& ant_sys_conf);

        inline std::array<float, NUM_ELEMENT_3D_GRID> compute_lcs2gcs(const ArrayOrientationConfig& ort_conf) {

            // Compute matrix that transforms vector from Local to 
            // Global Coordinate System (LCS -> GCS)
            // Based on 3GPP TR 38.901 Section 7.1.3 and TR 36.873 Section 5.1.3

            /*
              INPUT:
                orientation = [bearing (deg), downtilt (deg), slant (deg)]
            
                3GPP defines the composite rotation matrix as:
                  R = Rz(alpha) * Ry(beta) * Rx(gamma)
              where:
                alpha = bearing (rotation about Z-axis)
                beta  = downtilt (rotation about Y-axis)
                gamma = slant   (rotation about X-axis)
            
              ******************************!NOTE!******************************
              3GPP describes this R as a rotation that "maps a vector in the GCS into 
              the LCS". In other words, they interpret R as rotating the *coordinate 
              frame*.
            
                BUT: This implementation interprets R as an **active rotation** —
              that is, rotating the vector itself from LCS to its equivalent in GCS.
            
                This is mathematically valid because R is an orthogonal matrix:
                  R^-1 = R^T => inverse and transpose yield the same result.
            
                Therefore:
                - 3GPP’s "frame rotation" (passive) using R
                - is equivalent to this function’s "vector rotation" (active) using R
                - just applied in the opposite conceptual sense.
            */

            const float a = ort_conf.alpha;
            const float b = ort_conf.beta;
            const float g = ort_conf.gamma;

            const float ca = std::cos(DEG2RAD(a));
            const float cb = std::cos(DEG2RAD(b));
            const float cg = std::cos(DEG2RAD(g));

            const float sa = std::sin(DEG2RAD(a));
            const float sb = std::sin(DEG2RAD(b));
            const float sg = std::sin(DEG2RAD(g));

            std::array<float, NUM_ELEMENT_3D_GRID> r;

            // Column 0
            r[0] =  ca*cb;
            r[1] =  sa*cb;
            r[2] = -sb;

            // Column 1
            r[3] =  ca*sb*sg - sa*cg;
            r[4] =  sa*sb*sg + ca*cg;
            r[5] =  cb*sg;

            // Column 2
            r[6] =  ca*sb*cg + sa*sg;
            r[7] =  sa*sb*cg - ca*sg;
            r[8] =  cb*cg;

            return r;
        }

        inline std::array<float, NUM_ELEMENT_3D_GRID> multiply_3x3_matrices_BtAB(
            const float* A,
            const float* B
        ) {
            std::array<float, NUM_ELEMENT_3D_GRID> T;
            std::array<float, NUM_ELEMENT_3D_GRID> M;

            for (std::size_t j = 0; j < 3; ++j) {
                const std::size_t j3 = 3*j;
                const std::size_t j3_1 = j3 + 1;
                const std::size_t j3_2 = j3 + 2;
                T[j3] = A[0]*B[j3] + A[3]*B[j3_1] + A[6]*B[j3_2];
                T[j3_1] = A[1]*B[j3] + A[4]*B[j3_1] + A[7]*B[j3_2];
                T[j3_2] = A[2]*B[j3] + A[5]*B[j3_1] + A[8]*B[j3_2];
            }

            for (std::size_t j = 0; j < 3; ++j) {
                const std::size_t j3 = 3*j;
                const std::size_t j3_1 = j3 + 1;
                const std::size_t j3_2 = j3 + 2;
                M[j3] = B[0]*T[j3] + B[1]*T[j3_1] + B[2]*T[j3_2];
                M[j3_1] = B[3]*T[j3] + B[4]*T[j3_1] + B[5]*T[j3_2];
                M[j3_2] = B[6]*T[j3] + B[7]*T[j3_1] + B[8]*T[j3_2];
            }

            return M;
        }

        namespace los {

            inline SphericalDirection get_spherical_unit_vector(const float phi, const float theta) {

                // Get spherical unit vector only for LoS component

                const float theta_rad = DEG2RAD(theta);
                const float phi_rad = DEG2RAD(phi);

                SphericalDirection sph_dir_conf;

                sph_dir_conf.sinphi = std::sin(phi_rad);
                sph_dir_conf.cosphi = std::cos(phi_rad);
                sph_dir_conf.sintheta = std::sin(theta_rad);
                sph_dir_conf.costheta = std::cos(theta_rad);

                sph_dir_conf.rhat[0] = sph_dir_conf.sintheta*sph_dir_conf.cosphi;
                sph_dir_conf.rhat[1] = sph_dir_conf.sintheta*sph_dir_conf.sinphi;
                sph_dir_conf.rhat[2] = sph_dir_conf.costheta;
                return sph_dir_conf;
            }
        } // namespace los

        namespace nlos {

            inline SphericalDirections get_spherical_unit_vectors(const std::vector<float>& phi, const std::vector<float>& theta) {

                // Get spherical unit vector only for NLoS component

                SphericalDirections sph_dir_conf;

                sph_dir_conf.rhat.resize(NUM_3D_AXIS*phi.size());
                sph_dir_conf.sintheta.resize(phi.size());
                sph_dir_conf.costheta.resize(phi.size());
                sph_dir_conf.sinphi.resize(phi.size());
                sph_dir_conf.cosphi.resize(phi.size());

                float* rhat_ptr = sph_dir_conf.rhat.data();
                float* st_ptr = sph_dir_conf.sintheta.data();
                float* ct_ptr = sph_dir_conf.costheta.data();
                float* sp_ptr = sph_dir_conf.sinphi.data();
                float* cp_ptr = sph_dir_conf.cosphi.data();

                for (const auto& [p, t] : std::views::zip(phi, theta)) {

                    const float t_rad = DEG2RAD(t);
                    const float p_rad = DEG2RAD(p);

                    *st_ptr = std::sin(t_rad);;
                    *ct_ptr = std::cos(t_rad);
                    *sp_ptr = std::sin(p_rad);
                    *cp_ptr = std::cos(p_rad);

                    *rhat_ptr++ = (*st_ptr)*(*cp_ptr);
                    *rhat_ptr++ = (*st_ptr)*(*sp_ptr);
                    *rhat_ptr++ = *ct_ptr;
                    
                    ++st_ptr;
                    ++ct_ptr;
                    ++sp_ptr;
                    ++st_ptr;
                }
                return sph_dir_conf;
            }
        } // namespace nlos
    }

    namespace field_pattern {

        namespace los {
            void compute_field_term(
                float* field_term,
                const antenna::AntennaArrayConfig& ant_array_conf,
                const geometry::los::SphericalDirection sph_dir_conf,
                const std::size_t ant_idx
            );

            void get_location_term(
                std::complex<float>* loc_term,
                const std::array<float, geometry::NUM_3D_AXIS>& rhat,
                const std::vector<float>& dbar,
                const float lambda_0,
                const std::size_t ant_idx
            );
        } // namespace los

        namespace nlos {
            void compute_field_term(
                float* field_term,
                const antenna::AntennaArrayConfig& ant_array_conf,
                const float theta,
                const float phi,
                const std::size_t ant_idx
            );
        }

        float get_polarization_angle(
            std::size_t antennaIdx,
            const AntennaArrayConfig& ant_array_conf
        );

    }
}