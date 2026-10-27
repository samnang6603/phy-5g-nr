#include <vector>
#include <array>
#include <ranges>
#include <fftw3.h>
#include "generator.hpp"
#include "../cluster_profiles/profile.hpp"
#include "../../channel_filter/channel_filter.hpp"
#include "../../../modules/carrier/carrier.hpp"
#include "../../../common/channel_processing/channel_processing.hpp"

/***************** Constants ************************/
constexpr float OFFSET_IN_SYMBOL = 0.5f;
constexpr std::size_t NORMAL_CYLIC_PREFIX_LENGTH = 14;
constexpr std::array<std::size_t, NORMAL_CYLIC_PREFIX_LENGTH> NORMAL_CYCLIC_PREFIX_ARRAY = {
    80, 72, 72, 72, 72, 72, 72,
    80, 72, 72, 72, 72, 72 ,72
};

namespace channels::cdl::response {

    static void get_path_filters();

    void generate_ofdm_response(
        OFDMResponseConfig& ofdm_response_conf,
        const TimeVaryingPathGainConfig& time_varying_conf,
        const pdp::DelayProfileConfig& pdp_conf,
        const antenna::AntennaSystemConfig& ant_sys_conf,
        const ChannelControlConfig& channel_ctrl_conf
    ) {
        // Generate OFDM response.

        filter::ChannelFilterConfig channel_filt_conf;

        const std::size_t num_paths = pdp_conf.ClusterTypes.size();
        std::vector<float> path_delays(num_paths);
        const auto& pdp_table = pdp_conf.Table;

        for (std::size_t i = 0; i < num_paths; ++i) {
            path_delays[i] = pdp_table[i].normalized_delay;
        }

        //std::copy(pdp_table.begin()->normalized_delay,
        //          pdp_table.end()->normalized_delay,
        //        path_delays.begin()
        //        );
        
        filter::design_channel_filter(channel_filt_conf, path_delays);

        
        const int delay = common::channel_processing::estimate_delay(
            time_varying_conf.ClusterResponse, 
            channel_filt_conf.path_filter_coeffs, 
            time_varying_conf.TimeVector.size(), 
            pdp_conf.ClusterTypes.size(), 
            ant_sys_conf.NumInputSignals,
            ant_sys_conf.NumOutputSignals
        );

        const modules::Carrier& carrier_conf = ofdm_response_conf.Carrier;
        const std::size_t nfft = carrier_conf.NumFFT;
        const float fs = carrier_conf.SampleRate;
        const std::uint64_t num_time_samples = channel_ctrl_conf.NumTimeSamples;

        // Number of samples in each slot
        const float sample_offset_in_sym = nfft*OFFSET_IN_SYMBOL;
        std::vector<std::size_t> sym_len(NORMAL_CYLIC_PREFIX_LENGTH);
        std::size_t samples_per_slot = 0;
        for (const auto& [v, u] : std::views::zip(NORMAL_CYCLIC_PREFIX_ARRAY,sym_len)) {
            u = v + nfft;
            samples_per_slot += u;
        }

        // TODO LONG-TERM: If shift in slot #, implement here as circularly shift by slot #
        // std::size_t slot_idx

        // Calculate full subframe range with number of time samples
        const std::size_t samples_per_sf = fs/1000;
        const std::size_t num_full_sf = static_cast<std::size_t>(num_time_samples)/samples_per_sf;

        // Calculate full slot range with number of time samples
        const std::size_t num_slots = carrier_conf.SlotsPerFrame()*num_full_sf;
        const std::size_t rem_samples = num_time_samples % samples_per_sf;
        // TODO: when rem_samples > 0, fit it with another slot
        






    }

}