#include <cstddef>
#include <limits>
#include <stdexcept>
#include <span>
#include <vector>
#include <cmath>
#include "cdl_pdp.hpp"
#include "cdl_pdp_lut.hpp"
#include "../cdl.hpp"

/***************** Constants ************************/



/***************** Function Implementations ************************/

namespace channels::cdl::pdp {

    void initializeDelayProfile(DelayProfileConfig& pdp_conf) {

        // Initialize delay profile necessary properties

        std::span<const CDLCluster> source_pdp_table;
        bool hasLoS;
        const PDP profile = pdp_conf.DelayProfile;

        switch (profile) {
            case PDP::CDL_A:
                source_pdp_table = CDL_A_PDP;
                hasLoS = false;
                break;

            case PDP::CDL_B:
                source_pdp_table = CDL_B_PDP;
                hasLoS = false;
                break;

            case PDP::CDL_C:
                source_pdp_table = CDL_C_PDP;
                hasLoS = false;
                break;

            case PDP::CDL_D:
                source_pdp_table = CDL_D_PDP;
                hasLoS = true;
                break;

            case PDP::CDL_E:
                source_pdp_table = CDL_E_PDP;
                hasLoS = true;
                break;

            default:
                throw std::invalid_argument(
                    "initializeDelayProfile >>> Invalid/Unsupported CDL profile"
                );
        }

        // Scale KFactor if profile has LOS
        float desired_KFactor;
        if (hasLoS) {
            desired_KFactor = pdp_conf.KFactor.value_or(std::numeric_limits<float>::quiet_NaN());
        } else {
            desired_KFactor = std::numeric_limits<float>::quiet_NaN();
        }

        std::vector<CDLCluster> pdp_table(
            source_pdp_table.begin(),
            source_pdp_table.end()
        );







    }

    static void scale_Delays_KFactor(
        std::vector<CDLCluster>& pdp_table,
        const float delay_spread,
        const float desired_KFactor
    ) {
        // Scale delays and KFactor

        if (!isnan(desired_KFactor)) {
            // calculate K-Factor from pdp, lookup power index
            // see Equation 7.7.6-2

            float total_power_lin = 0.0f;

            for (std::size_t i = 1; i < pdp_table.size(); ++i) {
                const auto& power = pdp_table[i].power_db;
                total_power_lin += DB2POW(power);
            }

            float K_model = pdp_table[0].power_db - POW2DB(total_power_lin);

            // scale the power of all taps excluding the 1st tap (0-delay
            // component) by subtract the difference of target K-Factor and the
            // calculated K-Factor
            // See Equation 7.7.6-1
            for (std::size_t i = 1; i < pdp_table.size(); ++i) {
                auto& power = pdp_table[i].power_db;
                power += K_model - desired_KFactor;  
            }

            // calculate the RMS delay spread after the K-factor adjustment
            float sum_p = 0.0f;
            float w1 = 0.0f;
            float w2 = 0.0f;
            for (const auto& v : pdp_table) {
                const float p = DB2POW(v.power_db);
                const float t = v.normalized_delay;
                sum_p += p;
                float tmp0 = p*t;
                w1 += tmp0;
                w2 += tmp0*t;
            }
            w1 /= sum_p;
            w2 /= sum_p;
            const float tau_rms_inv = 1.0f/std::sqrt(w2 - w1*w1);

            for (auto& v : pdp_table) {
                v.normalized_delay *= tau_rms_inv; 
            }

        } else {

        }


    }

}
