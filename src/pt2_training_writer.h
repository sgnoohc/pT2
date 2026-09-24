#ifndef PT2_TRAINING_WRITER_H
#define PT2_TRAINING_WRITER_H

#include <array>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "histograms.h"
#include "pt2.h"
#include "rootReader.h"

// Flat tree of pT2 candidate features for NN training (read by pt2_ml/train_*.py): one row per pT2
class Pt2TrainingWriter
{
public:
    Pt2TrainingWriter(const std::string &path, const HistogramManager &hists);
    ~Pt2TrainingWriter();

    // Set the event index stored with the following rows
    void beginEvent(Long64_t ievt);
    // Write one row for a selected pT2
    void add(const rootReader &reader, const pT2 &pt2);
    // Write the tree and close the file
    void close();

    const std::string &path() const { return path_; }

private:
    std::string path_;
    TFile *file_ = nullptr;
    TTree *tree_ = nullptr;

    Long64_t event_idx_ = -1;
    float ls_pt_, ls_eta_, ls_phi_;
    float pls_pt_, pls_eta_, pls_phi_;
    float pls_charge_, pls_nhit_;
    float delta_pt_, delta_eta_, delta_phi_, delta_R_;
    float md0_dxy_, md0_dz_, md1_dxy_, md1_dz_;
    float md0_rz_, md1_rz_, log_abs_md0_dxy_;
    float lst_dPhi_, betaIn_, betaOut_, dBeta_;
    float lst_zResGeo_, lst_zResKin_, dAngle_;
    int is_real_;
    std::array<int, kNCat> is_cat_; // one-hot layer-connection category
};

#endif
