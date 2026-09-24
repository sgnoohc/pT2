#include "pt2_training_writer.h"

#include <cmath>
#include <stdexcept>

Pt2TrainingWriter::Pt2TrainingWriter(const std::string &path, const HistogramManager &hists)
    : path_(path)
{
    file_ = new TFile(path_.c_str(), "RECREATE");
    if (file_->IsZombie()) throw std::runtime_error("cannot create " + path_);
    tree_ = new TTree("tree", "pT2 Training Data");

    tree_->Branch("event_idx", &event_idx_);
    tree_->Branch("ls_pt", &ls_pt_);
    tree_->Branch("ls_eta", &ls_eta_);
    tree_->Branch("ls_phi", &ls_phi_);
    tree_->Branch("pls_pt", &pls_pt_);
    tree_->Branch("pls_eta", &pls_eta_);
    tree_->Branch("pls_phi", &pls_phi_);
    tree_->Branch("pls_charge", &pls_charge_);
    tree_->Branch("pls_nhit", &pls_nhit_);
    tree_->Branch("pt2_delta_pt", &delta_pt_);
    tree_->Branch("pt2_delta_eta", &delta_eta_);
    tree_->Branch("pt2_delta_phi", &delta_phi_);
    tree_->Branch("pt2_delta_R", &delta_R_);
    tree_->Branch("pt2_md0_dxy", &md0_dxy_);
    tree_->Branch("pt2_md0_dz", &md0_dz_);
    tree_->Branch("pt2_md1_dxy", &md1_dxy_);
    tree_->Branch("pt2_md1_dz", &md1_dz_);
    tree_->Branch("pt2_md0_rz", &md0_rz_);
    tree_->Branch("pt2_md1_rz", &md1_rz_);
    tree_->Branch("log_abs_md0_dxy", &log_abs_md0_dxy_);
    tree_->Branch("lst_dPhi", &lst_dPhi_);
    tree_->Branch("betaIn", &betaIn_);
    tree_->Branch("betaOut", &betaOut_);
    tree_->Branch("dBeta", &dBeta_);
    tree_->Branch("lst_zResGeo", &lst_zResGeo_);
    tree_->Branch("lst_zResKin", &lst_zResKin_);
    tree_->Branch("dAngle", &dAngle_);
    tree_->Branch("is_real", &is_real_);
    // Same order as extra_cuts::getConnectionIndex
    for (int ci = 0; ci < kNCat; ++ci)
        tree_->Branch(("is_" + hists.catNames[ci]).c_str(), &is_cat_[ci]);
}

Pt2TrainingWriter::~Pt2TrainingWriter()
{
    delete file_; // also deletes tree_
}

void Pt2TrainingWriter::beginEvent(Long64_t ievt)
{
    event_idx_ = ievt;
}

void Pt2TrainingWriter::add(const rootReader &reader, const pT2 &pt2)
{
    size_t p = pt2.pls_idx, l = pt2.ls_idx;

    // phi is stored raw; the training script takes sin/cos itself
    ls_pt_ = reader.ls_pt->at(l);
    ls_eta_ = reader.ls_eta->at(l);
    ls_phi_ = reader.ls_phi->at(l);
    pls_pt_ = reader.pls_pt->at(p);
    pls_eta_ = reader.pls_eta->at(p);
    pls_phi_ = reader.pls_phi->at(p);
    pls_charge_ = reader.pls_charge->at(p);
    pls_nhit_ = reader.pls_nhit->at(p);

    delta_pt_ = pt2.delta_pt;
    delta_eta_ = pt2.delta_eta;
    delta_phi_ = pt2.delta_phi;
    delta_R_ = pt2.delta_r;

    md0_dxy_ = pt2.heli[0];
    md0_dz_ = pt2.heli[1];
    md1_dxy_ = pt2.heli[2];
    md1_dz_ = pt2.heli[3];
    md0_rz_ = pt2.rz_simple.first;
    md1_rz_ = pt2.rz_simple.second;
    log_abs_md0_dxy_ = std::log(std::abs(static_cast<float>(pt2.heli[0])) + 1e-6f);

    lst_dPhi_ = pt2.lst_delta_phi;
    betaIn_ = pt2.beta_in;
    betaOut_ = pt2.beta_out;
    dBeta_ = pt2.delta_beta;
    lst_zResGeo_ = pt2.z_res_geo;
    lst_zResKin_ = pt2.z_res_kin;
    dAngle_ = pt2.delta_angle;

    is_real_ = pt2.is_real ? 1 : 0;
    for (int ci = 0; ci < kNCat; ++ci)
        is_cat_[ci] = (pt2.combo_idx == ci) ? 1 : 0;

    tree_->Fill();
}

void Pt2TrainingWriter::close()
{
    file_->cd();
    tree_->Write();
    file_->Close();
}
