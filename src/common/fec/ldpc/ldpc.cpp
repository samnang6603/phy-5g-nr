#include "ldpc.hpp"
#include "ldpc_lut.hpp"
#include <stdexcept>
#include <vector>

/********************** Constants ************************/
static constexpr uint64_t BG_SELECTION_A_THRESHOLD_1 = 292;
static constexpr uint64_t BG_SELECTION_A_THRESHOLD_2 = 3824;
static constexpr uint64_t BG2_MAX_A = 3840;
static constexpr uint64_t BG1_MAX_A = 8448;
static constexpr float BG_SELECTION_R_THRESHOLD_1 = 0.25f;
static constexpr float BG_SELECTION_R_THRESHOLD_2 = 0.67f;
static constexpr uint8_t BG_SHIFT_TABLE_NUM_COL = 10;
static constexpr uint8_t BG_LAYER_COL = 2;
static constexpr uint64_t ZC_MAX = 384;

/*********************** Alias ****************************/
using ShiftTableCol = std::array<uint16_t, BG_SHIFT_TABLE_NUM_COL>;

namespace common::fec::ldpc {

    /******************* Function Declarations **************************/
    static BGN select_base_graph(const uint64_t A, const float R);
    static uint8_t select_Kb(const BGN BGn, const uint64_t A);
    static void select_lift_size_set_idx(
        uint16_t& Zc, 
        uint8_t& set_idx, 
        const uint64_t A,
        const uint8_t Kb
    );
    static void setup_edges_and_shifts_lut(
        std::vector<uint16_t>& edges,
        std::vector<uint16_t>& shifts,
        const BGN BGn,
        const uint16_t Zc,
        const uint8_t set_idx
    );
    static void setup_layer(
        std::vector<uint16_t>& layers,
        const std::vector<uint16_t> edges,
        const BGN BGn
    );

    /********************* Class constructor ******************************/
    nrLDPC::nrLDPC(const uint64_t M, const float R) 
        : A_(static_cast<uint64_t>(M*R)) {

        BGn_ = select_base_graph(A_, R);
        const uint8_t kb = select_Kb(BGn_, A_);
        select_lift_size_set_idx(Zc_, set_idx_, A_, kb);

        switch (BGn_) {
            case BGN::BG1:
                K_ = Zc_*BG1_MSG_CHUNKS;
                N_ = Zc_*BG1_ALL_CHUNKS;
                break;
            
            case BGN::BG2:
                K_ = Zc_*BG2_MSG_CHUNKS;
                N_ = Zc_*BG2_ALL_CHUNKS;
                break;
            
            default:
                throw std::runtime_error(
                    "nrLDPC::nrLDPC >>> K_N_calculations : Invalid base graph"
                );
        }

        n_F_ = K_ - A_;
        setup_edges_and_shifts_lut(edges_, shifts_, BGn_, Zc_, set_idx_);
        setup_layer(layers_,edges_,BGn_);
    }

    /******************* Helper functions **************************/
    static BGN select_base_graph(const uint64_t A, const float R) {
        // Select base graph according to info length A and rate R 
        // from 3GPP 38.212 7.2.2 LDPC base graph selection

        if (
            A <= BG_SELECTION_A_THRESHOLD_1 ||
            (A <= BG_SELECTION_A_THRESHOLD_2 && R <= BG_SELECTION_R_THRESHOLD_2) ||
            R <= BG_SELECTION_R_THRESHOLD_1
        ) {
            if (A > BG2_MAX_A) {
                throw std::runtime_error(
                    "nrLDPC::select_base_graph >>> BG2 requires segmentation for A > 3840"
                );
            }
            return BGN::BG2;
        }

        if (A > BG1_MAX_A) {
            throw std::runtime_error(
                "nrLDPC::select_base_graph >>> BG1 requires segmentation for A > 8448"
            );
        }

        return BGN::BG1;
    }

    static uint8_t select_Kb(const BGN BGn, const uint64_t A) {
        // Select Kb from 3GPP 38.212 section 5.2.2

        if (BGn == BGN::BG1) {
            return 22;
        } else {
            if (A > 640) {
                return 10;
            } else if (A > 560) {
                return 9;
            } else if (A > 192) {
                return 8;
            } else {
                return 6;
            }
        }
        throw std::runtime_error("nrLDPC::select_Kb >>> Kb not selected");
    }

    static void select_lift_size_set_idx(
        uint16_t& Zc, 
        uint8_t& set_idx, 
        const uint64_t A,
        const uint8_t Kb
    ) {
        // Select lifting size Zc and set index from 3GPP Table 5.3.2-1

        // initiate a candidate Zc
        Zc = ZC_MAX;
        uint64_t Zc_tmp;
        
        for (std::size_t i = 0; i < LIFTING_SETS.size(); ++i) {

            for (std::size_t j = 0; j < LIFTING_SETS[i].size(); ++j) {

                Zc_tmp = LIFTING_SETS[i][j];
                const uint64_t tmp = static_cast<uint64_t>(Kb)*static_cast<uint64_t>(Zc_tmp);
                if (tmp == A) {
                    Zc = Zc_tmp;
                    set_idx = i;
                    return;
                } else if (tmp > A) {
                    if (Zc_tmp < Zc) {
                        Zc = Zc_tmp;
                        set_idx = i;
                    }
                }
            }
        }

        //throw std::runtime_error(
        //    "nrLDPC::select_lift_size_set_idx >> Zc not selected"
        //);
    }

    static void setup_edges_and_shifts_lut(
        std::vector<uint16_t>& edges,
        std::vector<uint16_t>& shifts,
        const BGN BGn,
        const uint16_t Zc,
        const uint8_t set_idx
    ) {
        // Select edges and shifts from 3GPP TS 38.212 Table 5.3.2-2/3: LDPC base graph
        // The first two columns give the row and column index of the coeffcients in the exponent parity matrx.
        // The 3rd-10th colums give the shift coefficeints to set index 0-7 respectively.
        // edges: shift connections from check nodes to variable nodes 
        // shifts: shifts for each correponding edge

        const std::size_t set_col = set_idx + 2; // total 3 columns
        const auto* table = static_cast<const ShiftTableCol*>(nullptr);
        std::size_t num_rows = 0;

        switch (BGn) {
            case BGN::BG1:
                table = BG1_SHIFT_TABLE.data();
                num_rows = BG1_SHIFT_TABLE.size();
                break;

            case BGN::BG2:
                table = BG2_SHIFT_TABLE.data();
                num_rows = BG2_SHIFT_TABLE.size();
                break;

            default:
                throw std::runtime_error(
                    "nrLDPC::select_edges_and_shifts_lut >>> Invalid base shifts and edges LUT"
                );
        }

        // The *2 because the edges house row and column index of the coeffcients in
        // the exponent parity matrx
        edges.resize(num_rows*2); 
        shifts.resize(num_rows);

        for (std::size_t i = 0; i < num_rows; ++i) {
            const auto& t = table[i];
            edges[i]            = t[0];
            edges[i + num_rows] = t[1];
            shifts[i]           = t[set_col] % Zc;
        }
    }

    static void setup_layer(
        std::vector<uint16_t>& layers,
        const std::vector<uint16_t> edges,
        const BGN BGn
    ) {
        // Creates layer where each layer is defined by edges specified 
        // by start edge index and number of edges

        const std::size_t num_edge_row = edges.size()/2;
        uint8_t num_layer_row;

        switch (BGn) {
            case BGN::BG1:
                num_layer_row = BG1_LAYER_ROWS;
                break;

            case BGN::BG2:
                num_layer_row = BG2_LAYER_ROWS;
                break;

            default:
                throw std::runtime_error(
                    "nrLDPC::setup_layer >>> Invalid base graph for layer setup"
                );
        }

        layers.resize(BG_LAYER_COL*num_layer_row);
        std::size_t jLayer = 0;
        std::size_t iEdge = 0;

        do {
            layers[jLayer] = iEdge;
            layers[jLayer + num_layer_row] = 0;

            while (edges[iEdge] == jLayer) {
                ++layers[jLayer + num_layer_row];
                ++iEdge;

                if (iEdge == num_edge_row) {
                    break;
                }
            }
            ++jLayer;
            
        } while (iEdge < num_edge_row);
    }
}