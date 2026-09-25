#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include "initial_phase.hpp"
#include "../antenna_array/antenna.hpp"

constexpr uint8_t NUM_ANGLE_POLARIZATION_COMBO = 4;
constexpr float TWO_PI = static_cast<float>(2.0f*M_PI);

namespace channels::cdl::phase {

    static float calculate_d3D(
        const antenna::geometry::PositionConfig& tx_pos_conf,
        const antenna::geometry::PositionConfig& rx_pos_conf
    );

    std::vector<float> generateInitialPhase(
        pdp::DelayProfileConfig& pdp_conf,
        RandomStreamConfig& randstream_conf,
        const ChannelControlConfig& control_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf
    ) {
        // Generates initial channel phase

        // Phi Size: The 4 refers to 4 polarization combinations
        // (theta-theta,theta-phi,phi-theta,phi-phi), theta: zenith, phi: azimuth
        // See TR 38.901 section 7.5 for more details
        std::size_t cluster_len = pdp_conf.Table.size();
        std::vector<float> phi(cluster_len*NUMBER_OF_RAYS_PER_CLUSTER*NUM_ANGLE_POLARIZATION_COMBO);

        auto& randomStream = randstream_conf.Stream;

        bool isRandomInitPhase = control_conf.InitialPhase == "Random";

        if (isRandomInitPhase) {
            for (auto& v : phi) {
                v = static_cast<float>(randomStream.rand())*TWO_PI - M_PI;
            }
        } else {
            // If initial phase is a scalar
        }

        if (pdp_conf.HasLoS) {
            // See TR 38.901 Equation 7.5-29
            // Phase of exponential term with d_3D
            const float lambda_0 = ant_sys_conf.lambda_0;

            // Calculate d_3D
            const float d_3D = calculate_d3D(
                ant_sys_conf.TransmitAntennaArray.Position,
                ant_sys_conf.ReceiveAntennaArray.Position
            );

            // Update LoS component
            phi[0] = -TWO_PI*d_3D/lambda_0;

        }

        return phi;
    }

    static float calculate_d3D(
        const antenna::geometry::PositionConfig& tx_pos_conf,
        const antenna::geometry::PositionConfig& rx_pos_conf
    ) {
        // Calculates the euclidean distance between Tx and Rx positions

        /*
            d_3D: effective projection of element separation in the ray's
            direction of arrival/departure
            Imagine the ray as a laser beam coming from a certain direction 
            (described by zenith and azimuth angles). 
            Then imagine your antenna element sitting somewhere in space.
            Now:
            Drop a perpendicular from that antenna element onto the laser 
            beam's direction vector. Measure that distance along the 
            ray direction. That's d_3D.
            It tells you:
            How far "along" the ray this element is located.
            Which in turn tells you how much phase shift this particular ray 
            will accumulate by the time it hits (or leaves) this element.
            Thus d_3D is a 3D Euclidean distance between position of Tx and
            Rx: d_3D = ||r_Rx - r_Tx|| norm-2
        */

        const float xtx = tx_pos_conf.x;
        const float ytx = tx_pos_conf.y;
        const float ztx = tx_pos_conf.z;

        const float xrx = rx_pos_conf.x;
        const float yrx = rx_pos_conf.y;
        const float zrx = rx_pos_conf.z;

        const float dx = xtx - xrx;
        const float dy = ytx - yrx;
        const float dz = ztx - zrx;

        const float dsum = dx + dy + dz;
        return std::sqrtf(dsum*dsum);
    }
}