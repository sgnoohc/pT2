#include "pt2_ntuple_writer.h"

#include <stdexcept>

Pt2NtupleWriter::Pt2NtupleWriter(const std::string& path, rootReader& reader)
    : path_(path), reader_(reader)
{
    file_ = new TFile(path_.c_str(), "RECREATE");
    if (file_->IsZombie()) throw std::runtime_error("cannot create " + path_);

    // Copy the input branch structure, but no entries
    tree_ = reader_.inputTree->CloneTree(0);

    tree_->Branch("pT2_pt", &pt_);
    tree_->Branch("pT2_eta", &eta_);
    tree_->Branch("pT2_phi", &phi_);
    tree_->Branch("pT2_plsIdx", &plsIdx_);
    tree_->Branch("pT2_lsIdx", &lsIdx_);
    tree_->Branch("pT2_isFake", &isFake_);
    tree_->Branch("pT2_isUsed", &isUsed_);
    tree_->Branch("pT2_isDuplicate", &isDuplicate_);
    tree_->Branch("pT2_deltaPt", &deltaPt_);
    tree_->Branch("pT2_deltaEta", &deltaEta_);
    tree_->Branch("pT2_deltaPhi", &deltaPhi_);
    tree_->Branch("pT2_dR", &dR_);
    tree_->Branch("pT2_NNscore", &nnScore_);
    tree_->Branch("pT2_matched_simIdx", &matchedSimIdx_);
    tree_->Branch("sim_pT2_matched", &simMatched_);
}

Pt2NtupleWriter::~Pt2NtupleWriter()
{
    delete file_;  // also deletes tree_
}

void Pt2NtupleWriter::beginEvent()
{
    pt_.clear();
    eta_.clear();
    phi_.clear();
    plsIdx_.clear();
    lsIdx_.clear();
    isFake_.clear();
    isUsed_.clear();
    isDuplicate_.clear();
    deltaPt_.clear();
    deltaEta_.clear();
    deltaPhi_.clear();
    dR_.clear();
    nnScore_.clear();
    matchedSimIdx_.clear();

    if (reader_.sim_pt) simMatched_.assign(reader_.sim_pt->size(), 0);
}

void Pt2NtupleWriter::add(const pT2& pt2)
{
    size_t plsIdx = pt2.pls_idx;

    // Approximate overall pT2 kinematics using its internal pLS
    pt_.push_back(reader_.pls_pt->at(plsIdx));
    eta_.push_back(reader_.pls_eta->at(plsIdx));
    phi_.push_back(reader_.pls_phi->at(plsIdx));

    plsIdx_.push_back(plsIdx);
    lsIdx_.push_back(pt2.ls_idx);
    isFake_.push_back(!pt2.is_real);
    isUsed_.push_back(pt2.is_used);

    deltaPt_.push_back(pt2.delta_pt);
    deltaEta_.push_back(pt2.delta_eta);
    deltaPhi_.push_back(pt2.delta_phi);
    dR_.push_back(pt2.delta_r);
    nnScore_.push_back(pt2.nn_score);

    if (pt2.is_real) {
        int simIdx = reader_.pls_simIdx->at(plsIdx);
        matchedSimIdx_.push_back({simIdx});
        if (simIdx >= 0 && simIdx < (int)simMatched_.size()) simMatched_[simIdx] += 1;
    } else {
        matchedSimIdx_.push_back({});
    }
}

void Pt2NtupleWriter::endEvent()
{
    flagDuplicates();
    injectIntoTCs();
    tree_->Fill();
}

void Pt2NtupleWriter::close()
{
    file_->cd();
    tree_->Write();
    file_->Close();
}

// A pT2 is a duplicate if its sim track was matched by more than one pT2
void Pt2NtupleWriter::flagDuplicates()
{
    for (const auto& simIdxs : matchedSimIdx_) {
        bool isDup = false;
        if (!simIdxs.empty()) {
            int simIdx = simIdxs[0];
            isDup = simIdx >= 0 && simIdx < (int)simMatched_.size() && simMatched_[simIdx] > 1;
        }
        isDuplicate_.push_back(isDup ? 1 : 0);
    }
}

void Pt2NtupleWriter::injectIntoTCs()
{
    rootReader& r = reader_;
    if (!r.tc_pt || !r.sim_tcIdx || !r.sim_tcIdxAll)
        throw std::runtime_error("tc_ or sim_tc branches are not loaded in rootReader");

    for (size_t i = 0; i < pt_.size(); i++) {
        int newTcIdx = r.tc_pt->size();
        int primarySimIdx = matchedSimIdx_[i].empty() ? -1 : matchedSimIdx_[i][0];

        r.tc_pt->push_back(pt_[i]);
        r.tc_eta->push_back(eta_[i]);
        r.tc_phi->push_back(phi_[i]);
        r.tc_isFake->push_back(isFake_[i]);
        r.tc_isDuplicate->push_back(isDuplicate_[i]);
        r.tc_type->push_back(2);  // Spoof as pLS to pass the filter
        r.tc_simIdx->push_back(primarySimIdx);
        r.tc_simIdxAll->push_back(matchedSimIdx_[i]);
        r.tc_nhits->push_back(2);
        r.tc_nlayers->push_back(2);

        // pT2 has no Outer Tracker hits
        if (r.tc_nhitOT) r.tc_nhitOT->push_back(0);
        // Real pT2s are fully matched, fakes not at all
        if (r.tc_pMatched) r.tc_pMatched->push_back(isFake_[i] ? 0.0 : 1.0);
        // Must be the same length as tc_simIdxAll
        if (r.tc_simIdxAllFrac) {
            std::vector<float> frac;
            if (!matchedSimIdx_[i].empty()) frac.push_back(1.0);
            r.tc_simIdxAllFrac->push_back(frac);
        }

        // pT2 has no higher-level LST objects, only the underlying pLS
        if (r.tc_pt5Idx) r.tc_pt5Idx->push_back(-1);
        if (r.tc_pt3Idx) r.tc_pt3Idx->push_back(-1);
        if (r.tc_t5Idx) r.tc_t5Idx->push_back(-1);
        if (r.tc_plsIdx) r.tc_plsIdx->push_back(plsIdx_[i]);

        // Reverse truth match: tell the sim particle this track found it
        if (primarySimIdx >= 0 && primarySimIdx < (int)r.sim_tcIdxAll->size()) {
            r.sim_tcIdxAll->at(primarySimIdx).push_back(newTcIdx);
            if (r.sim_tcIdx->at(primarySimIdx) == -1) r.sim_tcIdx->at(primarySimIdx) = newTcIdx;
        }
    }
}
