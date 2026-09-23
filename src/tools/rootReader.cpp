#include "rootReader.h"

rootReader::rootReader() : inputFile(nullptr), inputTree(nullptr) {}

rootReader::~rootReader() {
    if (inputFile) inputFile->Close();
}

bool rootReader::Init(const std::string& filename, const std::string& treename) {
    inputFile = TFile::Open(filename.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) return false;
    
    inputFile->GetObject(treename.c_str(), inputTree);
    if (!inputTree) return false;

    // LS Branches
    inputTree->SetBranchAddress("ls_pt", &ls_pt);
    inputTree->SetBranchAddress("ls_eta", &ls_eta);
    inputTree->SetBranchAddress("ls_phi", &ls_phi);
    inputTree->SetBranchAddress("ls_isFake", &ls_isFake);
    inputTree->SetBranchAddress("ls_simIdx", &ls_simIdx);
    inputTree->SetBranchAddress("ls_mdIdx0", &ls_mdIdx0);
    inputTree->SetBranchAddress("ls_mdIdx1", &ls_mdIdx1);
    inputTree->SetBranchAddress("ls_simIdxAllFrac", &ls_simIdxAllFrac);
    // Track Candidate Branches
    inputTree->SetBranchAddress("tc_pt5Idx", &tc_pt5Idx);
    inputTree->SetBranchAddress("tc_pt3Idx", &tc_pt3Idx);
    inputTree->SetBranchAddress("tc_t5Idx", &tc_t5Idx);
    inputTree->SetBranchAddress("tc_pt", &tc_pt);
    inputTree->SetBranchAddress("tc_eta", &tc_eta);
    inputTree->SetBranchAddress("tc_phi", &tc_phi);
    inputTree->SetBranchAddress("tc_isFake", &tc_isFake);
    inputTree->SetBranchAddress("tc_isDuplicate", &tc_isDuplicate);
    inputTree->SetBranchAddress("tc_type", &tc_type);
    inputTree->SetBranchAddress("tc_simIdxAll", &tc_simIdxAll);
    inputTree->SetBranchAddress("tc_simIdx", &tc_simIdx);
    inputTree->SetBranchAddress("tc_nhits", &tc_nhits);
    inputTree->SetBranchAddress("tc_nlayers", &tc_nlayers);
    inputTree->SetBranchAddress("tc_nhitOT", &tc_nhitOT);
    inputTree->SetBranchAddress("tc_pMatched", &tc_pMatched);
    inputTree->SetBranchAddress("tc_simIdxAllFrac", &tc_simIdxAllFrac);
    inputTree->SetBranchAddress("tc_plsIdx", &tc_plsIdx);
    // PT5 Branches
    inputTree->SetBranchAddress("pT5_t5Idx", &pt5_t5Idx);
    inputTree->SetBranchAddress("pT5_plsIdx", &pt5_plsIdx);
    // T5 Branches
    inputTree->SetBranchAddress("t5_t3Idx0", &t5_t3Idx0);
    inputTree->SetBranchAddress("t5_t3Idx1", &t5_t3Idx1);
    // T3 Branches
    inputTree->SetBranchAddress("t3_lsIdx0", &t3_lsIdx0);
    inputTree->SetBranchAddress("t3_lsIdx1", &t3_lsIdx1);
    // PT3 Branches
    inputTree->SetBranchAddress("pT3_t3Idx", &pt3_t3Idx);
    inputTree->SetBranchAddress("pT3_plsIdx", &pt3_plsIdx);
    // PLS Branches
    inputTree->SetBranchAddress("pLS_pt", &pls_pt);
    inputTree->SetBranchAddress("pLS_eta", &pls_eta);
    inputTree->SetBranchAddress("pLS_phi", &pls_phi);
    inputTree->SetBranchAddress("pLS_charge", &pls_charge);
    inputTree->SetBranchAddress("pLS_nhit", &pls_nhit);
    inputTree->SetBranchAddress("pLS_simIdxAllFrac", &pls_simIdxAllFrac);
    inputTree->SetBranchAddress("pLS_isFake", &pls_isFake);
    inputTree->SetBranchAddress("pLS_simIdx", &pls_simIdx);
    inputTree->SetBranchAddress("pLS_hit0_x", &pls_hit0_x);
    inputTree->SetBranchAddress("pLS_hit0_y", &pls_hit0_y);
    inputTree->SetBranchAddress("pLS_hit0_z", &pls_hit0_z);
    inputTree->SetBranchAddress("pLS_hit1_x", &pls_hit1_x);
    inputTree->SetBranchAddress("pLS_hit1_y", &pls_hit1_y);
    inputTree->SetBranchAddress("pLS_hit1_z", &pls_hit1_z);
    inputTree->SetBranchAddress("pLS_hit2_x", &pls_hit2_x);
    inputTree->SetBranchAddress("pLS_hit2_y", &pls_hit2_y);
    inputTree->SetBranchAddress("pLS_hit2_z", &pls_hit2_z);
    inputTree->SetBranchAddress("pLS_hit3_x", &pls_hit3_x);
    inputTree->SetBranchAddress("pLS_hit3_y", &pls_hit3_y);
    inputTree->SetBranchAddress("pLS_hit3_z", &pls_hit3_z);
    inputTree->SetBranchAddress("pLS_circleCenterX", &pLS_circleCenterX);
    inputTree->SetBranchAddress("pLS_circleCenterY", &pLS_circleCenterY);
    // Mini Doublet Branches
    inputTree->SetBranchAddress("md_anchor_x", &md_anchor_x);
    inputTree->SetBranchAddress("md_anchor_y", &md_anchor_y);
    inputTree->SetBranchAddress("md_anchor_z", &md_anchor_z);
    inputTree->SetBranchAddress("md_other_x", &md_other_x);
    inputTree->SetBranchAddress("md_other_y", &md_other_y);
    inputTree->SetBranchAddress("md_other_z", &md_other_z);
    inputTree->SetBranchAddress("md_layer", &md_layer);
    inputTree->SetBranchAddress("md_detId", &md_detId);
    // Sim Track Branches
    inputTree->SetBranchAddress("sim_pt", &sim_pt);
    inputTree->SetBranchAddress("sim_eta", &sim_eta);
    inputTree->SetBranchAddress("sim_phi", &sim_phi);
    inputTree->SetBranchAddress("sim_vx", &sim_vx);
    inputTree->SetBranchAddress("sim_vy", &sim_vy);
    inputTree->SetBranchAddress("sim_vz", &sim_vz);
    inputTree->SetBranchAddress("sim_pdgId", &sim_pdgId);
    inputTree->SetBranchAddress("sim_simHitX", &sim_simHitX);
    inputTree->SetBranchAddress("sim_simHitY", &sim_simHitY);
    inputTree->SetBranchAddress("sim_simHitZ", &sim_simHitZ);
    inputTree->SetBranchAddress("sim_tcIdx", &sim_tcIdx);
    inputTree->SetBranchAddress("sim_tcIdxAll", &sim_tcIdxAll);
    return true;
}

void rootReader::GetEntry(Long64_t i) {
    if (inputTree) inputTree->GetEntry(i);
}

Long64_t rootReader::GetEntries() const {
    return inputTree ? inputTree->GetEntries() : 0;
}


