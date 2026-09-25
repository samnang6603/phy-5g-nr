#include <cstddef>
#include <limits>
#include <stdexcept>
#include <span>
#include <vector>
#include <cmath>
#include "profile.hpp"
#include "predefined_profiles.hpp"
#include "../../../utils/utils.hpp"

/***************** Constants ************************/



/***************** Function Implementations ************************/

namespace channels::cdl::pdp {

    static void scale_KFactor(
        std::vector<CDLCluster>& pdp_table,
        const float desired_KFactor
    );

    static void scale_delay_spread(
        std::vector<CDLCluster>& pdp_table,
        const float delay_spread
    );
    
    static void get_per_cluster_param(
        DelayProfileConfig& pdp_conf,
        const CDLPerClusterParam per_cluster_param
    );

    static void get_cluster_types(
        DelayProfileConfig& pdp_conf,
        const std::span<const CDLCluster> source_pdp_table,
        const bool hasLoS
    );

    static void get_subcluster(
        DelayProfileConfig& pdp_conf
    );

    void initializeDelayProfile(DelayProfileConfig& pdp_conf) {

        // Initialize delay profile necessary properties

        std::span<const CDLCluster> source_pdp_table;
        std::vector<CDLCluster>& pdp_table = pdp_conf.Table;
        CDLPerClusterParam source_per_cluster_params;
        bool hasLoS;
        const float delay_spread = pdp_conf.DelaySpread;
        const PDP profile = pdp_conf.ProfileName;

        switch (profile) {
            case PDP::CDL_A:
                source_pdp_table = CDL_A_PDP;
                source_per_cluster_params = CDL_A_PER_CLUSTER_PARAM;
                hasLoS = false;
                break;

            case PDP::CDL_B:
                source_pdp_table = CDL_B_PDP;
                source_per_cluster_params = CDL_B_PER_CLUSTER_PARAM;
                hasLoS = false;
                break;

            case PDP::CDL_C:
                source_pdp_table = CDL_C_PDP;
                source_per_cluster_params = CDL_C_PER_CLUSTER_PARAM;
                hasLoS = false;
                break;

            case PDP::CDL_D:
                source_pdp_table = CDL_D_PDP;
                source_per_cluster_params = CDL_D_PER_CLUSTER_PARAM;
                hasLoS = true;
                break;

            case PDP::CDL_E:
                source_pdp_table = CDL_E_PDP;
                source_per_cluster_params = CDL_E_PER_CLUSTER_PARAM;
                hasLoS = true;
                break;

            default:
                throw std::invalid_argument(
                    "initializeDelayProfile >>> Invalid/Unsupported CDL profile"
                );
        }

        // Assign pdp table
        pdp_table.assign(
            source_pdp_table.begin(),
            source_pdp_table.end()
        );
        pdp_conf.HasLoS = hasLoS;

        // Build progation condition (cluster types) vector
        get_cluster_types(pdp_conf, source_pdp_table, hasLoS);

        // Get per cluster param
        get_per_cluster_param(pdp_conf, source_per_cluster_params);

        // Scale KFactor if profile has LOS
        if (hasLoS) {
            const float desired_KFactor = pdp_conf.KFactor.value_or(std::numeric_limits<float>::quiet_NaN());
            if (!isnan(desired_KFactor)) {
                scale_KFactor(pdp_table, desired_KFactor);
            }; 
        }

        // Scale dealy spread
        scale_delay_spread(pdp_table, delay_spread);


    }

    /*
    static void get_subcluster(
        DelayProfileConfig& pdp_conf
    ) {

        // TBI

    } */

    static void get_cluster_types(
        DelayProfileConfig& pdp_conf,
        const std::span<const CDLCluster> source_pdp_table,
        const bool hasLoS
    ) {
        // Get cluster types

        std::vector<PropagationCondition>& ClusterTypes = pdp_conf.ClusterTypes;

        ClusterTypes.assign(source_pdp_table.size(),PropagationCondition::NLOS);

        if (hasLoS) {
            ClusterTypes[0] = PropagationCondition::LOS;
        }
    }

    static void get_per_cluster_param(
        DelayProfileConfig& pdp_conf,
        const CDLPerClusterParam per_cluster_param
    ) {
        // Get corresponding pdp per cluster params: C_ASD, C_ASA, C_ZSD, C_ZSA, XPR

        pdp_conf.AngleSpreads.C_ASA = per_cluster_param.C_ASA;
        pdp_conf.AngleSpreads.C_ASD = per_cluster_param.C_ASD;
        pdp_conf.AngleSpreads.C_ZSA = per_cluster_param.C_ZSA;
        pdp_conf.AngleSpreads.C_ZSD = per_cluster_param.C_ZSD;
        pdp_conf.XPR = per_cluster_param.XPR;
    }

    static void scale_delay_spread(
        std::vector<CDLCluster>& pdp_table,
        const float delay_spread
    ) {
        // Scales pdp table by delay spread

        for (auto& v : pdp_table) {
            v.normalized_delay *= delay_spread;
        }
    }

    static void scale_KFactor(
        std::vector<CDLCluster>& pdp_table,
        const float desired_KFactor
    ) {
        // Scales KFactor
        // calculates K-Factor from pdp, lookup power index
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
    }
}
