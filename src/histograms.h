#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <TH1D.h>
#include <vector>
#include <string>
#include <utility>

constexpr int kNCat = 13;    // LS layer-connection categories
constexpr int kNCharge = 2;  // 0 = Pos, 1 = Neg

using HistGrid = TH1D* [kNCat][kNCharge];

// All per-pT2 variables for one (truth, usage, category, charge) slot
struct Pt2HistSet {
    TH1D *deltaPT, *deltaETA, *deltaPHI, *deltaR, *deltaAngle;
    TH1D *pls_ETA, *ls_ETA;
    TH1D *LSTdPhi, *LSTdBeta, *LSTbetaOut, *LSTOrgZRes, *LSTKinZRes;
    TH1D *MD0_dXY, *MD0_dZ, *MD1_dXY, *MD1_dZ;
    TH1D *MD0_rz_simple, *MD1_rz_simple;
};

class HistogramManager {
public:
    HistogramManager() = default;
    ~HistogramManager() = default;

    // Book empty histograms (detached from any TFile)
    void init();
    // Write all histograms to a new ROOT file
    void write(const std::string& path);
    // Read all histograms back from a file made by write()
    void load(const std::string& path);

    // Histograms to fill for a pT2 of given truth / usage in (cat, charge)
    const Pt2HistSet& set(bool real, bool unused, int cat, int charge) const {
        return sets_[real ? 0 : 1][unused ? 1 : 0][cat][charge];
    }

    // --- The 13 Valid Connections ---
    std::vector<std::string> catNames = {
        "L1F_to_L2F", "L1F_to_L2T", "L1T_to_L2F", "L1T_to_L2T", "L1T_to_E1PS",
        "L2F_to_L3F", "L2F_to_L3T", "L2T_to_L3F", "L2T_to_L3T", "L2T_to_E1PS",
        "E1PS_to_E2PS", "E1PS_to_E22S", "E2PS_to_E3PS"
    };
    std::vector<std::string> catTitles = {
        "Bar L1F -> Bar L2F", "Bar L1F -> Bar L2T", "Bar L1T -> Bar L2F", "Bar L1T -> Bar L2T", "Bar L1T -> Enc D1(PS)",
        "Bar L2F -> Bar L3F", "Bar L2F -> Bar L3T+", "Bar L2T -> Bar L3F", "Bar L2T -> Bar L3T", "Bar L2T -> Enc D1(PS)",
        "Enc D1(PS) -> Enc D2(PS)", "Enc D1(PS) -> Enc D2(2S)", "Enc D2(PS) -> Enc D3(PS)"
    };

    std::vector<std::string> chargeNames = {"Pos", "Neg"};
    std::vector<std::string> chargeTitles = {"Positive", "Negative"};

    // =========================================================================
    // GENERAL VARIABLES (Change all of these to [13])
    // =========================================================================
    TH1D* real_pt2_deltaPT[13][2];
    TH1D* real_pt2_deltaETA[13][2];
    TH1D* real_pt2_deltaPHI[13][2];
    TH1D* real_pt2_deltaR[13][2];
    TH1D* real_pt2_deltaAngle[13][2];

    TH1D* fake_pt2_deltaPT[13][2];
    TH1D* fake_pt2_deltaETA[13][2];
    TH1D* fake_pt2_deltaPHI[13][2];
    TH1D* fake_pt2_deltaR[13][2];
    TH1D* fake_pt2_deltaAngle[13][2];

    TH1D* real_unused_pt2_deltaPT[13][2];
    TH1D* real_unused_pt2_deltaETA[13][2];
    TH1D* real_unused_pt2_deltaPHI[13][2];
    TH1D* real_unused_pt2_deltaR[13][2];
    TH1D* real_unused_pt2_deltaAngle[13][2];

    TH1D* fake_unused_pt2_deltaPT[13][2];
    TH1D* fake_unused_pt2_deltaETA[13][2];
    TH1D* fake_unused_pt2_deltaPHI[13][2];
    TH1D* fake_unused_pt2_deltaR[13][2];
    TH1D* fake_unused_pt2_deltaAngle[13][2];

    // --- pLS Absolute Eta ---
    TH1D* real_pt2_pls_ETA[13][2];
    TH1D* fake_pt2_pls_ETA[13][2];
    TH1D* real_unused_pt2_pls_ETA[13][2];
    TH1D* fake_unused_pt2_pls_ETA[13][2];

    // --- LS Absolute Eta ---
    TH1D* real_pt2_ls_ETA[13][2];
    TH1D* fake_pt2_ls_ETA[13][2];
    TH1D* real_unused_pt2_ls_ETA[13][2];
    TH1D* fake_unused_pt2_ls_ETA[13][2];

    // LST specific variables
    TH1D* real_pt2_LSTdPhi[13][2];
    TH1D* fake_pt2_LSTdPhi[13][2];
    TH1D* real_unused_pt2_LSTdPhi[13][2];
    TH1D* fake_unused_pt2_LSTdPhi[13][2];

    TH1D* real_pt2_LSTdBeta[13][2];
    TH1D* fake_pt2_LSTdBeta[13][2];
    TH1D* real_unused_pt2_LSTdBeta[13][2];
    TH1D* fake_unused_pt2_LSTdBeta[13][2];

    TH1D* real_pt2_LSTbetaOut[13][2];
    TH1D* fake_pt2_LSTbetaOut[13][2];
    TH1D* real_unused_pt2_LSTbetaOut[13][2];
    TH1D* fake_unused_pt2_LSTbetaOut[13][2];

    TH1D* real_pt2_LSTKinZRes[13][2];
    TH1D* fake_pt2_LSTKinZRes[13][2];
    TH1D* real_unused_pt2_LSTKinZRes[13][2];
    TH1D* fake_unused_pt2_LSTKinZRes[13][2];

    TH1D* real_pt2_LSTOrgZRes[13][2];
    TH1D* fake_pt2_LSTOrgZRes[13][2];
    TH1D* real_unused_pt2_LSTOrgZRes[13][2];
    TH1D* fake_unused_pt2_LSTOrgZRes[13][2];

    // MD0 Components
    TH1D* real_pt2_MD0_dXY[13][2];
    TH1D* real_pt2_MD0_dZ[13][2];
    TH1D* fake_pt2_MD0_dXY[13][2];
    TH1D* fake_pt2_MD0_dZ[13][2];

    TH1D* real_unused_pt2_MD0_dXY[13][2];
    TH1D* real_unused_pt2_MD0_dZ[13][2];
    TH1D* fake_unused_pt2_MD0_dXY[13][2];
    TH1D* fake_unused_pt2_MD0_dZ[13][2];

    // MD1 Components
    TH1D* real_pt2_MD1_dXY[13][2];
    TH1D* real_pt2_MD1_dZ[13][2];
    TH1D* fake_pt2_MD1_dXY[13][2];
    TH1D* fake_pt2_MD1_dZ[13][2];

    TH1D* real_unused_pt2_MD1_dXY[13][2];
    TH1D* real_unused_pt2_MD1_dZ[13][2];
    TH1D* fake_unused_pt2_MD1_dXY[13][2];
    TH1D* fake_unused_pt2_MD1_dZ[13][2];

    // Separated R-Z Simple Pointing
    TH1D* real_pt2_MD0_rz_simple[13][2];
    TH1D* real_pt2_MD1_rz_simple[13][2];
    TH1D* fake_pt2_MD0_rz_simple[13][2];
    TH1D* fake_pt2_MD1_rz_simple[13][2];
    
    TH1D* real_unused_pt2_MD0_rz_simple[13][2];
    TH1D* real_unused_pt2_MD1_rz_simple[13][2];
    TH1D* fake_unused_pt2_MD0_rz_simple[13][2];
    TH1D* fake_unused_pt2_MD1_rz_simple[13][2];

private:
    // Every histogram grid above, keyed by its base name
    std::vector<std::pair<std::string, HistGrid*>> grids();
    void buildSets();

    Pt2HistSet sets_[2][2][kNCat][kNCharge];
};

#endif
