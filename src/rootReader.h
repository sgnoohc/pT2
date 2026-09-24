#ifndef ROOTREADER_H
#define ROOTREADER_H

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <iostream>

class rootReader {
public:
    rootReader();
    ~rootReader();
    bool Init(const std::string& filename, const std::string& treename);
    void GetEntry(Long64_t i);
    Long64_t GetEntries() const;

    // LS Variables
    std::vector<float> *ls_pt = nullptr;
    std::vector<float> *ls_eta = nullptr;
    std::vector<float> *ls_phi = nullptr;
    std::vector<std::vector<float>>* ls_simIdxAllFrac  = nullptr;
    std::vector<int>   *ls_isFake = nullptr;
    std::vector<int>   *ls_simIdx = nullptr;
    std::vector<int>   *ls_mdIdx0 = nullptr;
    std::vector<int>   *ls_mdIdx1 = nullptr;
    // TC Variables
    std::vector<int>   *tc_pt5Idx = nullptr;
    std::vector<int>   *tc_pt3Idx = nullptr;
    std::vector<int>   *tc_t5Idx = nullptr;
    std::vector<float> *tc_pt = nullptr;
    std::vector<float>* tc_eta = nullptr;
    std::vector<float>* tc_phi = nullptr;
    std::vector<int>* tc_isFake = nullptr;
    std::vector<int>* tc_isDuplicate = nullptr;
    std::vector<int>* tc_type = nullptr;
    std::vector<std::vector<int>>*  tc_simIdxAll = nullptr;
    std::vector<int>* tc_simIdx = nullptr;
    std::vector<int>* tc_nhits = nullptr;
    std::vector<int>* tc_nlayers = nullptr;
    std::vector<int>* tc_nhitOT = nullptr;
    std::vector<float>* tc_pMatched = nullptr;
    std::vector<std::vector<float>>* tc_simIdxAllFrac = nullptr;
    std::vector<int>* tc_plsIdx = nullptr;
    // pT5 Variables
    std::vector<int>   *pt5_plsIdx = nullptr;
    std::vector<int>   *pt5_t5Idx = nullptr;
    // T5 Variables
    std::vector<int>   *t5_t3Idx0 = nullptr;
    std::vector<int>   *t5_t3Idx1 = nullptr;
    // T3 Variables
    std::vector<int>   *t3_lsIdx0 = nullptr;
    std::vector<int>   *t3_lsIdx1 = nullptr;
    // pT3 Variables
    std::vector<int>   *pt3_t3Idx = nullptr;
    std::vector<int>   *pt3_plsIdx = nullptr;
    // pLS Variables
    std::vector<float> *pls_pt = nullptr;
    std::vector<float> *pls_eta = nullptr;
    std::vector<float> *pls_phi = nullptr;
    std::vector<int>   *pls_charge = nullptr;
    std::vector<int>   *pls_nhit = nullptr;
    std::vector<int>   *pls_isFake = nullptr;
    std::vector<int>   *pls_simIdx = nullptr; 
    std::vector<float> *pls_hit0_x = nullptr;
    std::vector<float> *pls_hit0_y = nullptr;
    std::vector<float> *pls_hit0_z = nullptr;
    std::vector<float> *pls_hit1_x = nullptr;
    std::vector<float> *pls_hit1_y = nullptr;
    std::vector<float> *pls_hit1_z = nullptr;
    std::vector<float> *pls_hit2_x = nullptr;
    std::vector<float> *pls_hit2_y = nullptr;
    std::vector<float> *pls_hit2_z = nullptr;
    std::vector<float> *pls_hit3_x = nullptr;
    std::vector<float> *pls_hit3_y = nullptr;
    std::vector<float> *pls_hit3_z = nullptr;
    std::vector<std::vector<float>>* pls_simIdxAllFrac = nullptr;
    std::vector<float>* pLS_circleCenterX = nullptr;
    std::vector<float>* pLS_circleCenterY = nullptr;
    // Mini-Doublet Variables
    std::vector<int>   *md_layer = nullptr;
    std::vector<int>   *md_detId = nullptr;
    std::vector<float> *md_anchor_x = nullptr;
    std::vector<float> *md_anchor_y = nullptr;
    std::vector<float> *md_anchor_z = nullptr;
    std::vector<float> *md_other_x = nullptr;
    std::vector<float> *md_other_y = nullptr;
    std::vector<float> *md_other_z = nullptr;
    // Sim Variables
    std::vector<float> *sim_pt = nullptr;    
    std::vector<float>* sim_eta = nullptr;
    std::vector<float>* sim_phi = nullptr;
    std::vector<float>* sim_vx = nullptr;
    std::vector<float>* sim_vy = nullptr;
    std::vector<float>* sim_vz = nullptr;
    std::vector<int>* sim_pdgId = nullptr;
    std::vector<int>* sim_tcIdx = nullptr;
    std::vector<std::vector<int>>* sim_tcIdxAll = nullptr;
    std::vector<std::vector<float>>* sim_simHitX = nullptr;
    std::vector<std::vector<float>>* sim_simHitY = nullptr;
    std::vector<std::vector<float>>* sim_simHitZ = nullptr;

    // Computed Variables
    std::vector<float> pls_origin_z;
    std::vector<int> pls_superbin;
    std::vector<bool> pls_isUsed;
    std::vector<bool> ls_isUsed;

private:
    TFile* inputFile;

public:
    TTree* inputTree;
};
#endif
