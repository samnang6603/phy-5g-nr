#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include "cdl_phase.hpp"
#include "../cdl.hpp"
#include "../../matlab_mt19937.hpp"

constexpr uint8_t NUM_ANGLE_POLARIZATION_COMBO = 4;
constexpr float TWO_PI = static_cast<float>(2.0f*M_PI);

namespace channels::cdl::phase {

    std::vector<float> generateInitialPhase(
        pdp::DelayProfileConfig& pdp_conf,
        nrCDLChannel::RandomStreamConfig& rstream_conf,
        nrCDLChannel::ChannelControlConfig& control_conf,
        antenna::AntennaSystemConfig& ant_conf
    ) {
        // Generates initial channel phase

        // Phi Size: The 4 refers to 4 polarization combinations
        // (theta-theta,theta-phi,phi-theta,phi-phi), theta: zenith, phi: azimuth
        // See TR 38.901 section 7.5 for more details
        std::size_t cluster_len = pdp_conf.Table.size();
        std::vector<float> phi(cluster_len*NUMBER_OF_RAYS*NUM_ANGLE_POLARIZATION_COMBO);

        std::uint32_t seed = rstream_conf.Seed;

        MATLABMT19937 randomStream(seed);

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
            const float lambda_0 = SPEED_OF_LIGHT/ant_conf.CarrierFrequency;

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

            //const float tx_pos_sum = ant_conf.TransmitAntennaArray;

            //const float d_3D = 1;

        }

        return phi;
    }
}