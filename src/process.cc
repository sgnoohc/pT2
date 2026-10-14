#include "process.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "extra_cuts.h"
#include "extrapolation.h"
#include "gator.h"
#include "histograms.h"
#include "pt2.h"
#include "pt2_cuts.h"
#include "pt2_ntuple_writer.h"
#include "pt2_scorer.h"
#include "pt2_training_writer.h"
#include "rootReader.h"
#include "tools.h"

namespace
{

    // Only pT2s whose pLS is in this pT window are studied
    constexpr float kPlsPtMin = 0.6;
    constexpr float kPlsPtMax = 0.8;

    // pT2s are scored by the NN in chunks of this size
    constexpr size_t kNNBatchSize = 4096;

    struct SuperbinMaps
    {
        SuperbinToDetIdMap pos, neg, non;
    };

    // Build all pT2 candidates of the current event
    void buildEventPt2s(rootReader &reader, const SuperbinMaps &maps, bool lowPT, pT2Collection &pt2s)
    {
        reader.pls_origin_z.clear();
        reader.pls_superbin.clear();
        pt2s.clear();

        UsedMask usedMask = buildUsedMask(reader);
        reader.ls_isUsed = std::move(usedMask.ls_isUsed);
        reader.pls_isUsed = std::move(usedMask.pls_isUsed);

        size_t nLS = reader.ls_pt->size();
        size_t nPLS = reader.pls_pt->size();

        DetIdToLSMap detidToLS;
        for (size_t k = 0; k < nLS; ++k)
            for (int detId : getDetIdsForLS(reader, k))
                detidToLS[detId].push_back(k);

        reader.pls_origin_z.reserve(nPLS);
        reader.pls_superbin.reserve(nPLS);
        for (size_t j = 0; j < nPLS; ++j)
        {
            reader.pls_origin_z.push_back(CalculatePlsZ(reader, j));
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j, lowPT));
            buildPt2sForPLS(j, reader, maps.pos, maps.neg, maps.non, detidToLS, pt2s);
        }
    }

    // Fill the pT2's truth / kinematic fields and features.
    // nn_score is left at -1; scoreBatch() sets it
    void computeFeatures(const rootReader &reader, pT2 &pt2)
    {
        size_t p = pt2.pls_idx, l = pt2.ls_idx;

        pt2.delta_pt = deltaPt(reader.pls_pt->at(p), reader.ls_pt->at(l));
        pt2.delta_eta = deltaEta(reader.pls_eta->at(p), reader.ls_eta->at(l));
        pt2.delta_phi = deltaPhi(reader.pls_phi->at(p), reader.ls_phi->at(l));
        pt2.is_real = pt2TruthFinder(reader, p, l);
        pt2.is_used = pt2UsedCalculator(reader, p, l);

        uint32_t detId0 = reader.md_detId->at(reader.ls_mdIdx0->at(l));
        uint32_t detId1 = reader.md_detId->at(reader.ls_mdIdx1->at(l));
        pt2.combo_idx = extra_cuts::getConnectionIndex(extra_cuts::getCategoryFromDetId(detId0), extra_cuts::getCategoryFromDetId(detId1));
        pt2.charge_idx = (reader.pls_charge->at(p) > 0) ? 0 : 1;

        pt2.delta_r = std::sqrt(pt2.delta_eta * pt2.delta_eta + pt2.delta_phi * pt2.delta_phi);
        pt2.pls_eta = reader.pls_eta->at(p);
        pt2.ls_eta = reader.ls_eta->at(l);

        pt2.heli = extrapolation::extrapolatePlsHelicallyAndGetDistance(p, l, reader);
        pt2.rz_simple = extrapolation::extrapolateSimplePointingInRZ(p, l, reader);
        pt2.delta_angle = extrapolation::calculateDeltaAngle(p, l, reader);

        pt2.lst_delta_phi = extra_cuts::calculateLSTDPhi(p, l, reader);
        std::vector<double> betas = extra_cuts::calculateLSTdBeta(p, l, reader);
        pt2.beta_in = betas[0];
        pt2.beta_out = betas[1];
        pt2.delta_beta = betas[2];
        pt2.z_res_geo = extra_cuts::calculateLSTOriginZResidual(p, l, reader);
        pt2.z_res_kin = extra_cuts::calculateLSTKinematicZResidual(p, l, reader);
    }

    // Set nn_score for a batch of pT2s with one NN call
    void scoreBatch(const Pt2Scorer &scorer, const rootReader &reader, const std::vector<pT2 *> &batch)
    {
        std::vector<Pt2Scorer::FeatureVector> inputs;
        inputs.reserve(batch.size());
        for (const pT2 *pt2 : batch) inputs.push_back(Pt2Scorer::features(reader, *pt2));

        std::vector<float> scores = scorer.score(inputs);
        for (size_t i = 0; i < batch.size(); ++i) batch[i]->nn_score = scores[i];
    }

    // Values at or below the sentinels mean "could not be computed" and are skipped
    void fillSet(const Pt2HistSet &h, const pT2 &pt2)
    {
        h.deltaPT->Fill(pt2.delta_pt);
        h.deltaETA->Fill(pt2.delta_eta);
        h.deltaPHI->Fill(pt2.delta_phi);
        h.deltaR->Fill(pt2.delta_r);
        h.pls_ETA->Fill(pt2.pls_eta);
        h.ls_ETA->Fill(pt2.ls_eta);

        if (pt2.delta_angle > -1.0) h.deltaAngle->Fill(pt2.delta_angle);
        if (pt2.lst_delta_phi > -100.0) h.LSTdPhi->Fill(pt2.lst_delta_phi);
        if (pt2.delta_beta > -100.0) h.LSTdBeta->Fill(pt2.delta_beta);
        if (pt2.beta_out > -100.0) h.LSTbetaOut->Fill(pt2.beta_out);
        if (pt2.z_res_geo > -100.0) h.LSTOrgZRes->Fill(pt2.z_res_geo);
        if (pt2.z_res_kin > -100.0) h.LSTKinZRes->Fill(pt2.z_res_kin);

        if (pt2.heli[0] >= 0)
        {
            h.MD0_dXY->Fill(pt2.heli[0]);
            h.MD0_dZ->Fill(pt2.heli[1]);
            h.MD1_dXY->Fill(pt2.heli[2]);
            h.MD1_dZ->Fill(pt2.heli[3]);
        }

        if (pt2.rz_simple.first > -900)
        {
            h.MD0_rz_simple->Fill(pt2.rz_simple.first);
            h.MD1_rz_simple->Fill(pt2.rz_simple.second);
        }
    }

    // Ambiguity resolution: several pT2s in the same zone often share a pLS and
    // are really the same track reconstructed more than once. Keep only the
    // `keep` best of each (pLS, zone, charge) group, ranked by how well the pLS
    // extrapolation lands on the LS, and drop the rest from `selected`.
    //
    // This runs after the threshold cuts, so a pT2 that was already rejected
    // never takes a slot away from one that survived.
    void resolveAmbiguity(std::vector<pT2 *> &selected, int keep)
    {
        std::map<long long, std::vector<std::pair<double, pT2 *>>> groups;
        for (pT2 *p : selected)
        {
            if (p->combo_idx < 0 || p->charge_idx < 0) continue;
            long long key = (long long)p->pls_idx * 26 + p->combo_idx * 2 + p->charge_idx;
            groups[key].push_back({pt2_cuts::ambiguityChi2(*p), p});
        }

        std::vector<pT2 *> losers;
        for (auto &kv : groups)
        {
            auto &v = kv.second;
            if ((int)v.size() <= keep) continue;
            std::nth_element(v.begin(), v.begin() + keep, v.end());
            for (size_t k = keep; k < v.size(); ++k)
                losers.push_back(v[k].second);
        }
        if (losers.empty()) return;

        std::sort(losers.begin(), losers.end());
        selected.erase(std::remove_if(selected.begin(), selected.end(),
                                      [&](pT2 *p)
                                      { return std::binary_search(losers.begin(), losers.end(), p); }),
                       selected.end());
    }

    void fillHistograms(const HistogramManager &hists, const pT2 &pt2)
    {
        if (pt2.combo_idx < 0) return;
        fillSet(hists.set(pt2.is_real, false, pt2.combo_idx, pt2.charge_idx), pt2);
        if (!pt2.is_used) fillSet(hists.set(pt2.is_real, true, pt2.combo_idx, pt2.charge_idx), pt2);
    }

} // namespace

int runProcess(const Config &cfg)
{
    std::filesystem::create_directories(cfg.outputDir);

    SuperbinMaps maps;
    loadSuperbinDetIdMap(cfg.pixelMapDir, maps.pos, maps.neg, maps.non);

    rootReader reader;
    if (!reader.Init(cfg.inputFile, "tree")) throw std::runtime_error("could not open ROOT file or TTree: " + cfg.inputFile);

    HistogramManager hists;
    hists.init();

    std::unique_ptr<Pt2NtupleWriter> writer;
    std::unique_ptr<Pt2TrainingWriter> trainWriter;
    std::unique_ptr<Pt2Scorer> scorer;
    if (cfg.writeRoot)
    {
        if (cfg.useNN)
            scorer = std::make_unique<Pt2Scorer>(cfg.nnModelDir + "/model.onnx", cfg.nnModelDir + "/mean.npy", cfg.nnModelDir + "/std.npy");
        writer = std::make_unique<Pt2NtupleWriter>(cfg.outputDir + "/LSTNtuple_with_pT2.root", reader);
        trainWriter = std::make_unique<Pt2TrainingWriter>(cfg.outputDir + "/pt2_training_data.root", hists);
    }

    Long64_t totalEntries = reader.GetEntries();
    Long64_t nEntries = (cfg.nEvents > 0 && cfg.nEvents < totalEntries) ? cfg.nEvents : totalEntries;

    print_creature();

    pT2Collection pt2s;
    for (Long64_t ievt = 0; ievt < nEntries; ++ievt)
    {
        reader.GetEntry(ievt);
        if (ievt % 2 == 0) printProgressBar(ievt, nEntries);

        if (writer) writer->beginEvent();
        if (trainWriter) trainWriter->beginEvent(ievt);

        buildEventPt2s(reader, maps, cfg.lowPT, pt2s);

        // Everything done with a pT2 once its features and NN score are known
        auto processPt2 = [&](const pT2 &pt2)
        {
            // pt2.print();

            // The per-zone threshold cuts (-c) and the ambiguity step (-a N)
            // have already been applied by the time a pT2 gets here; see the
            // selection loop below.
            //
            // NN cut, score thresholds by real efficiency:
            //   90%: 0.99761337  95%: 0.99330878  96%: 0.98942512  97%: 0.97682154
            //   98%: 0.92524022  99%: 0.71041596  99.9%: 0.01175442
            // if (pt2.nn_score < 0.92524022f) return;

            if (writer) writer->add(pt2);
            if (trainWriter) trainWriter->add(reader, pt2);
            fillHistograms(hists, pt2);
        };

        size_t npt2 = 0;   // pT2s of this event inside the pT window

        // Select the pT2s of this event: pT window, then the threshold cuts,
        // then ambiguity. All three are off by default, so without -c / -a this
        // keeps exactly what it always did.
        std::vector<pT2 *> selected;
        for (auto &pt2 : pt2s)
        {
            float plsPt = reader.pls_pt->at(pt2.pls_idx);
            if (!(plsPt > kPlsPtMin && plsPt < kPlsPtMax)) continue;

            npt2++;
            computeFeatures(reader, pt2);
            if (cfg.applyCuts && !pt2_cuts::passAll(pt2, reader.ls_pt->at(pt2.ls_idx))) continue;
            selected.push_back(&pt2);
        }
        if (cfg.keepPerPls > 0) resolveAmbiguity(selected, cfg.keepPerPls);

        // Score the survivors in batches, then process them
        std::vector<pT2 *> batch;
        auto flush = [&]()
        {
            if (batch.empty()) return;
            if (scorer) scoreBatch(*scorer, reader, batch);
            for (pT2 *pt2 : batch) processPt2(*pt2);
            batch.clear();
        };

        for (pT2 *pt2 : selected)
        {
            batch.push_back(pt2);
            if (batch.size() == kNNBatchSize) flush();
        }
        flush();

        std::cout <<  " npt2: " << npt2 <<  std::endl;

        if (writer) writer->endEvent();
    }
    std::cout << "\n";

    hists.write(cfg.histFile);
    std::cout << "Saved histograms to: " << cfg.histFile << std::endl;

    if (writer)
    {
        writer->close();
        std::cout << "Saved ROOT file containing old + new pT2 branches to: " << writer->path() << std::endl;
    }
    if (trainWriter)
    {
        trainWriter->close();
        std::cout << "Saved training data to: " << trainWriter->path() << std::endl;
    }

    return 0;
}
