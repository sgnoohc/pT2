#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <getopt.h>
#include <cmath>
#include <iomanip>

#include "gator.h"
#include "histograms.h"
#include "plotting.h"
#include "plot_recipes.h"
#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h"
#include "extra_cuts.h"
#include "model_inference.h"

#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv) {

    // Default Arguments
    bool makePlots = false;
    bool writeRoot = false;
    bool lowPT = false;
    int nEvents = -1;
    double targetPercent = 90;
    double myCutZ0 = 0.4896;  
    double myCutZ1 = 0.9304;
    std::string inputFile; 
    std::string outputDir;
    float minPt = 0.6;
    float maxPt = 0.8;

    // Command Line Arguments
    int opt;
    while ((opt = getopt(argc, argv, "prki:o:n:e:m:M:")) != -1) {
        switch (opt) {
            case 'p':
                makePlots = true;
                break;
            case 'r':
                writeRoot = true;
                break;
            case 'k':
                lowPT = true;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputDir = optarg;
                break;
            case 'n':
                nEvents = std::stoi(optarg);
                break; 
            case 'e':
                targetPercent = std::stod(optarg);
                break;
            case 'm':
                minPt = std::stof(optarg);
                break;
            case 'M':
                maxPt = std::stof(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0] << "\n"
                    << "[-p] Make Plots\n"
                    << "[-k] Run Low pT\n"
                    << "[-r] Make Root File \n"
                    << "[-i] Input File Path \n"
                    << "[-o] Output Directory \n"
                    << "[-n] Number of Events \n"
                    << "[-e] Target percent \n"
                    << "[-m GeV] Minimum pT\n"
                    << "[-M GeV] Maximum pT\n";
           return 1;
        }    
    }

    // Code should do something
    if ((!makePlots) && (!writeRoot)){
        std::cerr << "Error: the code should have some sort of output!" << std::endl;
        return 1;
    }
    // Set input and output file paths
    if (outputDir.empty()){
        outputDir = "output";
    }
    if (inputFile.empty()){ 
        if (lowPT){
            //inputFile = "/blue/avery/aaponteutani/CMSSW_16_1_0/src/RecoTracker/LSTCore/standalone/LSTNtuple_pionGun_0p3.root"; //0.5GeV
           // inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple_LowPT.root"; //0.6 GeV
           //   inputFile = "/blue/avery/aaponteutani/CMSSW_16_1_0/src/RecoTracker/LSTCore/standalone/LSTNtuple_PU200_0p5.root";//0.5 GeV
            inputFile = "/blue/avery/aaponteutani/CMSSW_16_1_0/src/RecoTracker/LSTCore/standalone/LSTNtuple_PU200_0p6.root";//0.6 GeV
        }
        else {
            inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple.root";
        }
    } 


    // Print Gator
    print_gator();

    // Print Configuration
    std::cout << "\n=== Configuration ===\n";
    std::cout << "Input file:        " << inputFile << "\n";
    std::cout << "Output directory:  " << outputDir << "\n";
    std::cout << "Make plots:        " << (makePlots ? "yes" : "no") << "\n";
    std::cout << "Use Low pT:        " << (lowPT ? "yes" : "no") << "\n";
    std::cout << "Write ROOT file:   " << (writeRoot ? "yes" : "no") << "\n";
    if (nEvents > 0) std::cout << "Number of events:  " << nEvents << "\n";
    std::cout << "=====================\n\n";

    // Create the output directory
    std::filesystem::create_directories(outputDir);

    const std::string MODEL_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/model.onnx";
    const std::string MEAN_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/mean.npy";
    const std::string STD_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/std.npy";
    // ----- CUT SCORES -----
    //const float CUT_SCORE = 0.9976133704185486f; // 90% Real Efficiency
    //const float CUT_SCORE = 0.9933087825775146f; // 95% Real Efficiency
    //const float CUT_SCORE = 0.9894251227378845f; // 96% Real Efficiency
    //const float CUT_SCORE = 0.9768215417861938f; // 97% Real Efficiency
    const float CUT_SCORE = 0.9252402186393738f; // 98% Real Efficiency
    //const float CUT_SCORE = 0.7104159593582153f; // 99% Real Efficiency
    //const float CUT_SCORE = 0.011754416860640049f; // 99.9% Real Efficiency

    Pt2Scorer scorer(MODEL_PATH, MEAN_PATH, STD_PATH);
    std::cout << "ML model loaded.  Active cut score = " << CUT_SCORE << "\n\n";

    // Initialize Histograms
    HistogramManager hists;
    hists.init();

    // Get the Correct Pixel Map Directory
    std::string pixelMapFileDir;
    if (lowPT){
       pixelMapFileDir = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Pixel_Maps_0p6GeV/";
       // pixelMapFileDir = "/blue/p.chang/aaponteutani/LSTGeometry/output_0p4/pixelmap/"; //0.4 GeV
       // pixelMapFileDir = "/blue/p.chang/aaponteutani/LSTGeometry/output_0p3/pixelmap/"; //0.5 GeV
    }
    else{
        pixelMapFileDir = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Pixel_Maps_0p8GeV/";
    }

    // Load superbin --> detID
    SuperbinToDetIdMap superbinToDetIds_POS;
    SuperbinToDetIdMap superbinToDetIds_NEG;
    SuperbinToDetIdMap superbinToDetIds_NON;
    try {
        loadSuperbinDetIdMap(pixelMapFileDir, 
                            superbinToDetIds_POS, 
                            superbinToDetIds_NEG,
                            superbinToDetIds_NON);
    }
    catch (const std::exception& e) {
        std::cerr << "Error loading superbin map: " << e.what() << std::endl;
        return 1;
    }

    // Initialize per-event variables
    DetIdToLSMap detidToLS;
    pT2Collection pt2s;

    // Setup the Root Reader
    rootReader reader;
    if (!reader.Init(inputFile, "tree")) {
        std::cerr << "Error: Could not find Root File or TTree! \n" << std::endl;
        return 1;
    } 

    Long64_t totalEntries = reader.GetEntries();
    Long64_t entriesToProcess = (nEvents > 0 && nEvents < totalEntries) ? nEvents : totalEntries;
    // Setting up training tree
    std::string trainFileName = outputDir + "/pt2_training_data.root";
    TFile* trainFile = new TFile(trainFileName.c_str(), "RECREATE");
    TTree* trainTree = new TTree("tree", "pT2 Training Data");

    // === Training Tree: per-event accumulator vectors ===
    const std::vector<std::string> kCatNames = {
        "E1PS_to_E22S",   // 0  — Endcap D1 PS  → Endcap D2 2S
        "E1PS_to_E2PS",   // 1  — Endcap D1 PS  → Endcap D2 PS
        "E2PS_to_E3PS",   // 2  — Endcap D2 PS  → Endcap D3 PS
        "L1F_to_L2F",     // 3  — Layer 1 Flat  → Layer 2 Flat
        "L1F_to_L2T",     // 4  — Layer 1 Flat  → Layer 2 Tilted
        "L1T_to_E1PS",    // 5  — Layer 1 Tilt  → Endcap D1 PS
        "L1T_to_L2F",     // 6  — Layer 1 Tilt  → Layer 2 Flat
        "L1T_to_L2T",     // 7  — Layer 1 Tilt  → Layer 2 Tilted
        "L2F_to_L3F",     // 8  — Layer 2 Flat  → Layer 3 Flat
        "L2F_to_L3T",     // 9  — Layer 2 Flat  → Layer 3 Tilted
        "L2T_to_E1PS",    // 10 — Layer 2 Tilt  → Endcap D1 PS
        "L2T_to_L3F",     // 11 — Layer 2 Tilt  → Layer 3 Flat
        "L2T_to_L3T",     // 12 — Layer 2 Tilt  → Layer 3 Tilted
    };
    std::vector<float> t_ls_pt, t_ls_eta, t_ls_phi;
    std::vector<float> t_pls_pt, t_pls_eta, t_pls_phi;
    std::vector<float> t_pls_charge, t_pls_nhit;
    std::vector<float> t_pt2_delta_pt, t_pt2_delta_eta, t_pt2_delta_phi, t_pt2_delta_R;
    std::vector<float> t_pt2_md0_dxy, t_pt2_md0_dz, t_pt2_md1_dxy, t_pt2_md1_dz;
    std::vector<float> t_pt2_md0_rz, t_pt2_md1_rz, t_log_abs_md0_dxy;
    std::vector<float> t_lst_dPhi, t_betaIn, t_betaOut, t_dBeta;
    std::vector<float> t_lst_zResGeo, t_lst_zResKin, t_dAngle;
    std::vector<int>   t_is_real;
    std::array<std::vector<int>, 13> t_is_cat; // one-module/type: is_E1PS_to_E22S, …

    trainTree->Branch("ls_pt",           &t_ls_pt);
    trainTree->Branch("ls_eta",          &t_ls_eta);
    trainTree->Branch("ls_phi",      &t_ls_phi);
    trainTree->Branch("pls_pt",          &t_pls_pt);
    trainTree->Branch("pls_eta",         &t_pls_eta);
    trainTree->Branch("pls_phi",     &t_pls_phi);
    trainTree->Branch("pls_charge",      &t_pls_charge);
    trainTree->Branch("pls_nhit",        &t_pls_nhit);
    trainTree->Branch("pt2_delta_pt",    &t_pt2_delta_pt);
    trainTree->Branch("pt2_delta_eta",   &t_pt2_delta_eta);
    trainTree->Branch("pt2_delta_phi",   &t_pt2_delta_phi);
    trainTree->Branch("pt2_delta_R",     &t_pt2_delta_R);
    trainTree->Branch("pt2_md0_dxy",     &t_pt2_md0_dxy);
    trainTree->Branch("pt2_md0_dz",      &t_pt2_md0_dz);
    trainTree->Branch("pt2_md1_dxy",     &t_pt2_md1_dxy);
    trainTree->Branch("pt2_md1_dz",      &t_pt2_md1_dz);
    trainTree->Branch("pt2_md0_rz",      &t_pt2_md0_rz);
    trainTree->Branch("pt2_md1_rz",      &t_pt2_md1_rz);
    trainTree->Branch("log_abs_md0_dxy", &t_log_abs_md0_dxy);
    trainTree->Branch("lst_dPhi",        &t_lst_dPhi);
    trainTree->Branch("betaIn",          &t_betaIn);
    trainTree->Branch("betaOut",         &t_betaOut);
    trainTree->Branch("dBeta",           &t_dBeta);
    trainTree->Branch("lst_zResGeo",     &t_lst_zResGeo);
    trainTree->Branch("lst_zResKin",     &t_lst_zResKin);
    trainTree->Branch("dAngle",          &t_dAngle);
    trainTree->Branch("is_real",         &t_is_real);
    for (int ci = 0; ci < 13; ++ci)
        trainTree->Branch(("is_" + kCatNames[ci]).c_str(), &t_is_cat[ci]);

    
    // === NEW: Output ROOT File Setup ===
    TFile* outRootFile = nullptr;
    TTree* outTree = nullptr;

    // Vectors to hold our new pT2 Branches
    std::vector<float> pT2_pt, pT2_eta, pT2_phi;
    std::vector<int>   pT2_plsIdx, pT2_lsIdx;
    std::vector<int>   pT2_isFake, pT2_isUsed, pT2_isDuplicate;
    std::vector<float> pT2_delta_pt, pT2_delta_eta, pT2_delta_phi;
    std::vector<float> pT2_dR;
    std::vector<std::vector<int>> pT2_matched_simIdx; 
    std::vector<int>   sim_pT2_matched;
    std::vector<float>   pT2_NNscore;
    // You can also add vector<vector<int>> for pT2_simIdxAll here if you calculate it!

    if (writeRoot) {
        std::string outRootName = outputDir + "/LSTNtuple_with_pT2.root";
        outRootFile = new TFile(outRootName.c_str(), "RECREATE");

        // Clone the input tree structure. '0' means copy branches but not entries yet.
        // NOTE: Make sure your rootReader makes the TTree pointer accessible.
        // Often it's named 'tree' or retrieved via a getter like 'GetTree()'.
        outTree = reader.inputTree->CloneTree(0);

        // Attach our new branches
        outTree->Branch("pT2_pt",      &pT2_pt);
        outTree->Branch("pT2_eta",     &pT2_eta);
        outTree->Branch("pT2_phi",     &pT2_phi);
        outTree->Branch("pT2_plsIdx",  &pT2_plsIdx);
        outTree->Branch("pT2_lsIdx",   &pT2_lsIdx);
        outTree->Branch("pT2_isFake",  &pT2_isFake);
        outTree->Branch("pT2_isUsed",  &pT2_isUsed);
        outTree->Branch("pT2_isDuplicate", &pT2_isDuplicate);
        outTree->Branch("pT2_deltaPt", &pT2_delta_pt);
        outTree->Branch("pT2_deltaEta",&pT2_delta_eta);
        outTree->Branch("pT2_deltaPhi",&pT2_delta_phi);
        outTree->Branch("pT2_dR",      &pT2_dR);
        outTree->Branch("pT2_matched_simIdx", &pT2_matched_simIdx); 
        outTree->Branch("sim_pT2_matched", &sim_pT2_matched);
        outTree->Branch("pT2_NNscore", &pT2_NNscore);
    }
    // ===================================

    print_creature();

    // Main Looper
    for (Long64_t ievt = 0; ievt < entriesToProcess; ++ievt) {

        reader.GetEntry(ievt);

        // Update progress bar every N events
        if (ievt % 2 == 0 || ievt == entriesToProcess)
            printProgressBar(ievt, entriesToProcess);

        // Clear Variables from the last event
        reader.pls_origin_z.clear();
        reader.pls_superbin.clear();
        detidToLS.clear();
        pt2s.clear();
        reader.pls_isUsed.clear();
        reader.ls_isUsed.clear();
        
        // === NEW: Clear pT2 vectors from previous event ===
        if (writeRoot) {
            pT2_pt.clear();  pT2_eta.clear(); pT2_phi.clear();
            pT2_plsIdx.clear(); pT2_lsIdx.clear();
            pT2_isFake.clear(); pT2_isUsed.clear(); pT2_isDuplicate.clear();
            pT2_delta_pt.clear(); pT2_delta_eta.clear(); pT2_delta_phi.clear();
            pT2_dR.clear();
            pT2_matched_simIdx.clear();
            pT2_NNscore.clear();

            // Clear training vectors
            t_ls_pt.clear(); t_ls_eta.clear(); t_ls_phi.clear(); 
            t_pls_pt.clear(); t_pls_eta.clear(); t_pls_phi.clear(); 
            t_pls_charge.clear(); t_pls_nhit.clear();
            t_pt2_delta_pt.clear(); t_pt2_delta_eta.clear(); t_pt2_delta_phi.clear(); t_pt2_delta_R.clear();
            t_pt2_md0_dxy.clear(); t_pt2_md0_dz.clear(); t_pt2_md1_dxy.clear(); t_pt2_md1_dz.clear();
            t_pt2_md0_rz.clear(); t_pt2_md1_rz.clear(); t_log_abs_md0_dxy.clear();
            t_lst_dPhi.clear(); t_betaIn.clear(); t_betaOut.clear(); t_dBeta.clear();
            t_lst_zResGeo.clear(); t_lst_zResKin.clear(); t_dAngle.clear();
            t_is_real.clear();
            for (auto& v : t_is_cat) v.clear();
            
            if (reader.sim_pt) {
                sim_pT2_matched.assign(reader.sim_pt->size(), 0);
            }
        }
        // ==================================================

        // Get used masks for the LS and pLS
        UsedMask usedMask = buildUsedMask(reader);
        reader.ls_isUsed  = std::move(usedMask.ls_isUsed);
        reader.pls_isUsed = std::move(usedMask.pls_isUsed);

        // Get Size of pLS and LS
        size_t nLS = reader.ls_pt->size();
        size_t nPLS = reader.pls_pt->size();

        // Reserve New Calculated Variables
        reader.pls_origin_z.reserve(nPLS);
        reader.pls_superbin.reserve(nPLS);

        // LS Loop
        for (size_t k = 0; k < nLS; ++k) {
            std::vector<int> detIds = getDetIdsForLS(reader, k);
                for (int detId : detIds) {
                    detidToLS[detId].push_back(k);
                }
        }

        // pLS Loop
        for (size_t j = 0; j < nPLS; ++j) {
            reader.pls_origin_z.push_back(CalculatePlsZ(reader, j));
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j, lowPT));
            buildPt2sForPLS(j, reader, superbinToDetIds_POS, superbinToDetIds_NEG, superbinToDetIds_NON, detidToLS, pt2s);
        }

        // pT2 Loop
        for (auto& pt2 : pt2s) {
            size_t plsIdx = pt2.pls_idx;
            size_t lsIdx = pt2.ls_idx;
           
            // --- CATEGORIZE MD0 and MD1 ---
            int md0_idx = reader.ls_mdIdx0->at(lsIdx);
            int md1_idx = reader.ls_mdIdx1->at(lsIdx);
            uint32_t detId0 = reader.md_detId->at(md0_idx);
            uint32_t detId1 = reader.md_detId->at(md1_idx);
            
            int cat0 = extra_cuts::getCategoryFromDetId(detId0);
            int cat1 = extra_cuts::getCategoryFromDetId(detId1);
            int comboIdx = extra_cuts::getConnectionIndex(cat0, cat1); 
            //-----Charge-------
            int charge = reader.pls_charge->at(plsIdx);
            int cIdx = (charge > 0) ? 0 : 1;

            pt2.delta_pt  = deltaPt(reader.pls_pt->at(plsIdx), reader.ls_pt->at(lsIdx));
            pt2.delta_eta = deltaEta(reader.pls_eta->at(plsIdx), reader.ls_eta->at(lsIdx));
            pt2.delta_phi = deltaPhi(reader.pls_phi->at(plsIdx), reader.ls_phi->at(lsIdx));
            pt2.is_real = pt2TruthFinder(reader, plsIdx, lsIdx);
            pt2.is_used = pt2UsedCalculator(reader, plsIdx, lsIdx);

            // Using only pls_pt from 0.6 to 0.8
            float pt = reader.pls_pt->at(plsIdx);

            if (pt < minPt || pt > maxPt) {
                continue;
                }

            // Physics Calculations
            float dR = std::sqrt(pt2.delta_eta * pt2.delta_eta + pt2.delta_phi * pt2.delta_phi);
            std::vector<double> heli = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);
            std::pair<double, double> rz_simple = extrapolation::extrapolateSimplePointingInRZ(plsIdx, lsIdx, reader);
            double dAngle = extrapolation::calculateDeltaAngle(plsIdx, lsIdx, reader);
            
            double lst_dPhi = extra_cuts::calculateLSTDPhi(plsIdx, lsIdx, reader);
            std::vector<double> betas = extra_cuts::calculateLSTdBeta(plsIdx, lsIdx, reader);
            double dBeta = betas[2];
            double betaOut = betas[1];
            double betaIn = betas[0];
            double lst_zResGeo = extra_cuts::calculateLSTOriginZResidual(plsIdx, lsIdx, reader);
            double lst_zResKin = extra_cuts::calculateLSTKinematicZResidual(plsIdx, lsIdx, reader);

            // Preparing the 21 features for the NN scorer
            std::array<float, 21> features = {{
                reader.ls_pt->at(lsIdx),
                reader.ls_eta->at(lsIdx),
                std::sin(reader.ls_phi->at(lsIdx)),
                std::cos(reader.ls_phi->at(lsIdx)),
                reader.pls_pt->at(plsIdx),
                reader.pls_eta->at(plsIdx),
                std::sin(reader.pls_phi->at(plsIdx)),
                std::cos(reader.pls_phi->at(plsIdx)),
                static_cast<float>(reader.pls_charge->at(plsIdx)),
                static_cast<float>(reader.pls_nhit->at(plsIdx)),
                pt2.delta_pt,
                pt2.delta_eta,
                pt2.delta_phi,
                dR,
                static_cast<float>(heli[0]),
                static_cast<float>(heli[1]),
                static_cast<float>(heli[2]),
                static_cast<float>(heli[3]),
                static_cast<float>(rz_simple.first),
                static_cast<float>(rz_simple.second),
                std::log(std::abs(static_cast<float>(heli[0])) + 1e-6f),
            }};

            float mlscore = scorer.score(features);
            //Apply ML cut, uncomment if wanting to apply cuts!

            // if (
            //     mlscore >= CUT_SCORE &&
            //     true
            // )
            
            // === NEW: Push data to branches ===
            // This happens ONLY for pt2s passing the > 0.8 pT requirement above.
            if (writeRoot) {
                // Approximate overall pT2 kinematics using its internal pLS
                pT2_pt.push_back(reader.pls_pt->at(plsIdx));
                pT2_eta.push_back(reader.pls_eta->at(plsIdx));
                pT2_phi.push_back(reader.pls_phi->at(plsIdx));

                pT2_plsIdx.push_back(plsIdx);
                pT2_lsIdx.push_back(lsIdx);
                pT2_isFake.push_back(!pt2.is_real);
                pT2_isUsed.push_back(pt2.is_used);

                pT2_delta_pt.push_back(pt2.delta_pt);
                pT2_delta_eta.push_back(pt2.delta_eta);
                pT2_delta_phi.push_back(pt2.delta_phi);
                pT2_dR.push_back(dR);
                pT2_NNscore.push_back(mlscore);
                // -------- TRUTH MATCHING --------
                if (pt2.is_real) {
                    // Get the sim track ID from the pLS
                    int simIdx = reader.pls_simIdx->at(plsIdx);
                    pT2_matched_simIdx.push_back({simIdx});

                    // Register that this sim track was successfully found!
                    if (simIdx >= 0 && simIdx < sim_pT2_matched.size()) {
                        sim_pT2_matched[simIdx] += 1;
                    }
                } else {
                    // Push an empty vector for fakes
                    pT2_matched_simIdx.push_back({});
                }
                // === Fill training vectors (one entry per pT2 candidate) ===
                t_ls_pt.push_back(reader.ls_pt->at(lsIdx));
                t_ls_eta.push_back(reader.ls_eta->at(lsIdx));
                t_ls_phi.push_back(std::sin(reader.ls_phi->at(lsIdx)));
                t_pls_pt.push_back(reader.pls_pt->at(plsIdx));
                t_pls_eta.push_back(reader.pls_eta->at(plsIdx));
                t_pls_phi.push_back(std::sin(reader.pls_phi->at(plsIdx)));
                t_pls_charge.push_back(static_cast<float>(reader.pls_charge->at(plsIdx)));
                t_pls_nhit.push_back(static_cast<float>(reader.pls_nhit->at(plsIdx)));
                t_pt2_delta_pt.push_back(pt2.delta_pt);
                t_pt2_delta_eta.push_back(pt2.delta_eta);
                t_pt2_delta_phi.push_back(pt2.delta_phi);
                t_pt2_delta_R.push_back(dR);
                t_pt2_md0_dxy.push_back(static_cast<float>(heli[0]));
                t_pt2_md0_dz.push_back(static_cast<float>(heli[1]));
                t_pt2_md1_dxy.push_back(static_cast<float>(heli[2]));
                t_pt2_md1_dz.push_back(static_cast<float>(heli[3]));
                t_pt2_md0_rz.push_back(static_cast<float>(rz_simple.first));
                t_pt2_md1_rz.push_back(static_cast<float>(rz_simple.second));
                t_log_abs_md0_dxy.push_back(std::log(std::abs(static_cast<float>(heli[0])) + 1e-6f));
                t_lst_dPhi.push_back(static_cast<float>(lst_dPhi));
                t_betaIn.push_back(static_cast<float>(betaIn));
                t_betaOut.push_back(static_cast<float>(betaOut));
                t_dBeta.push_back(static_cast<float>(dBeta));
                t_lst_zResGeo.push_back(static_cast<float>(lst_zResGeo));
                t_lst_zResKin.push_back(static_cast<float>(lst_zResKin));
                t_dAngle.push_back(static_cast<float>(dAngle));
                t_is_real.push_back(pt2.is_real ? 1 : 0);
                // 1 for this pT2's category, 0 for all others
                for (int ci = 0; ci < 13; ++ci)
                    t_is_cat[ci].push_back((comboIdx == ci) ? 1 : 0);

            }

            //if(heli[1] > myCutZ0 ||  heli[3] > myCutZ1 ){continue;}
            //if(heli[0] > 2.3896 ||  heli[2] > 3.4234 ){continue;}
            //if(pt2.delta_phi < -0.3493 || pt2.delta_phi > 0.3457){continue;}
            //if(rz_simple.first < -2.8875  || rz_simple.first > 1.2586 || rz_simple.second < -4.7368 || rz_simple.second > 1.8688){continue;}
            //if(pt2.delta_pt < -0.6123 || pt2.delta_pt > 0.1846){continue;}
            //if(dBeta < -0.0445 || dBeta > 0.0393){continue;}
            //if(lst_zResKin < -3.5895 || lst_zResKin > 3.7145){continue;}
            if (comboIdx >= 0) {
            if (pt2.is_real) {
                hists.real_pt2_deltaPT[comboIdx][cIdx]->Fill(pt2.delta_pt);
                hists.real_pt2_deltaETA[comboIdx][cIdx]->Fill(pt2.delta_eta);
                hists.real_pt2_deltaPHI[comboIdx][cIdx]->Fill(pt2.delta_phi);
                hists.real_pt2_deltaR[comboIdx][cIdx]->Fill(dR);
                hists.real_pt2_pls_ETA[comboIdx][cIdx]->Fill(reader.pls_eta->at(plsIdx));
                hists.real_pt2_ls_ETA[comboIdx][cIdx]->Fill(reader.ls_eta->at(lsIdx));

                if (dAngle > -1.0) hists.real_pt2_deltaAngle[comboIdx][cIdx]->Fill(dAngle);
                if (lst_dPhi > -100.0) hists.real_pt2_LSTdPhi[comboIdx][cIdx]->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.real_pt2_LSTdBeta[comboIdx][cIdx]->Fill(dBeta);
                if (betaOut > -100.0) hists.real_pt2_LSTbetaOut[comboIdx][cIdx]->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.real_pt2_LSTOrgZRes[comboIdx][cIdx]->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.real_pt2_LSTKinZRes[comboIdx][cIdx]->Fill(lst_zResKin);

                // Fill Separated 3D components for Real
                if (heli[0] >= 0) {
                    hists.real_pt2_MD0_dXY[comboIdx][cIdx]->Fill(heli[0]); 
                    hists.real_pt2_MD0_dZ[comboIdx][cIdx]->Fill(heli[1]);
                    hists.real_pt2_MD1_dXY[comboIdx][cIdx]->Fill(heli[2]); 
                    hists.real_pt2_MD1_dZ[comboIdx][cIdx]->Fill(heli[3]);

                }

                if (rz_simple.first > -900) {
                    hists.real_pt2_MD0_rz_simple[comboIdx][cIdx]->Fill(rz_simple.first);
                    hists.real_pt2_MD1_rz_simple[comboIdx][cIdx]->Fill(rz_simple.second);
                }

                if (!pt2.is_used) {
                    hists.real_unused_pt2_deltaPT[comboIdx][cIdx]->Fill(pt2.delta_pt);
                    hists.real_unused_pt2_deltaETA[comboIdx][cIdx]->Fill(pt2.delta_eta);
                    hists.real_unused_pt2_deltaPHI[comboIdx][cIdx]->Fill(pt2.delta_phi);
                    hists.real_unused_pt2_deltaR[comboIdx][cIdx]->Fill(dR);
                    hists.real_unused_pt2_pls_ETA[comboIdx][cIdx]->Fill(reader.pls_eta->at(plsIdx));
                    hists.real_unused_pt2_ls_ETA[comboIdx][cIdx]->Fill(reader.ls_eta->at(lsIdx));

                    if (dAngle > -1.0) hists.real_unused_pt2_deltaAngle[comboIdx][cIdx]->Fill(dAngle);
                    if (lst_dPhi > -100.0) hists.real_unused_pt2_LSTdPhi[comboIdx][cIdx]->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.real_unused_pt2_LSTdBeta[comboIdx][cIdx]->Fill(dBeta);
                    if (betaOut > -100.0) hists.real_unused_pt2_LSTbetaOut[comboIdx][cIdx]->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.real_unused_pt2_LSTOrgZRes[comboIdx][cIdx]->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.real_unused_pt2_LSTKinZRes[comboIdx][cIdx]->Fill(lst_zResKin);

                    if (heli[0] >= 0) {
                        hists.real_unused_pt2_MD0_dXY[comboIdx][cIdx]->Fill(heli[0]); 
                        hists.real_unused_pt2_MD0_dZ[comboIdx][cIdx]->Fill(heli[1]);
                        hists.real_unused_pt2_MD1_dXY[comboIdx][cIdx]->Fill(heli[2]); 
                        hists.real_unused_pt2_MD1_dZ[comboIdx][cIdx]->Fill(heli[3]);
                    }

                    if (rz_simple.first > -900) {
                        hists.real_unused_pt2_MD0_rz_simple[comboIdx][cIdx]->Fill(rz_simple.first);
                        hists.real_unused_pt2_MD1_rz_simple[comboIdx][cIdx]->Fill(rz_simple.second);
                    }
                }
            } 
            else {
                hists.fake_pt2_deltaPT[comboIdx][cIdx]->Fill(pt2.delta_pt);
                hists.fake_pt2_deltaETA[comboIdx][cIdx]->Fill(pt2.delta_eta);
                hists.fake_pt2_deltaPHI[comboIdx][cIdx]->Fill(pt2.delta_phi);
                hists.fake_pt2_deltaR[comboIdx][cIdx]->Fill(dR);
                hists.fake_pt2_pls_ETA[comboIdx][cIdx]->Fill(reader.pls_eta->at(plsIdx));
                hists.fake_pt2_ls_ETA[comboIdx][cIdx]->Fill(reader.ls_eta->at(lsIdx));

                if (dAngle > -1.0) hists.fake_pt2_deltaAngle[comboIdx][cIdx]->Fill(dAngle);
                if (lst_dPhi > -100.0) hists.fake_pt2_LSTdPhi[comboIdx][cIdx]->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.fake_pt2_LSTdBeta[comboIdx][cIdx]->Fill(dBeta);
                if (betaOut > -100.0) hists.fake_pt2_LSTbetaOut[comboIdx][cIdx]->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.fake_pt2_LSTOrgZRes[comboIdx][cIdx]->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.fake_pt2_LSTKinZRes[comboIdx][cIdx]->Fill(lst_zResKin);

                // Fill Separated 3D components for Fake
                if (heli[0] >= 0) {
                    hists.fake_pt2_MD0_dXY[comboIdx][cIdx]->Fill(heli[0]); 
                    hists.fake_pt2_MD0_dZ[comboIdx][cIdx]->Fill(heli[1]);
                    hists.fake_pt2_MD1_dXY[comboIdx][cIdx]->Fill(heli[2]); 
                    hists.fake_pt2_MD1_dZ[comboIdx][cIdx]->Fill(heli[3]);

                }

                if (rz_simple.first > -900) {
                    hists.fake_pt2_MD0_rz_simple[comboIdx][cIdx]->Fill(rz_simple.first);
                    hists.fake_pt2_MD1_rz_simple[comboIdx][cIdx]->Fill(rz_simple.second);
                }

                if (!pt2.is_used) {
                    hists.fake_unused_pt2_deltaPT[comboIdx][cIdx]->Fill(pt2.delta_pt);
                    hists.fake_unused_pt2_deltaETA[comboIdx][cIdx]->Fill(pt2.delta_eta);
                    hists.fake_unused_pt2_deltaPHI[comboIdx][cIdx]->Fill(pt2.delta_phi);
                    hists.fake_unused_pt2_deltaR[comboIdx][cIdx]->Fill(dR);
                    hists.fake_unused_pt2_pls_ETA[comboIdx][cIdx]->Fill(reader.pls_eta->at(plsIdx));
                    hists.fake_unused_pt2_ls_ETA[comboIdx][cIdx]->Fill(reader.ls_eta->at(lsIdx));

                    if (dAngle > -1.0) hists.fake_unused_pt2_deltaAngle[comboIdx][cIdx]->Fill(dAngle);
                    if (lst_dPhi > -100.0) hists.fake_unused_pt2_LSTdPhi[comboIdx][cIdx]->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.fake_unused_pt2_LSTdBeta[comboIdx][cIdx]->Fill(dBeta);
                    if (betaOut > -100.0) hists.fake_unused_pt2_LSTbetaOut[comboIdx][cIdx]->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.fake_unused_pt2_LSTOrgZRes[comboIdx][cIdx]->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.fake_unused_pt2_LSTKinZRes[comboIdx][cIdx]->Fill(lst_zResKin);

                    if (heli[0] >= 0) {
                        hists.fake_unused_pt2_MD0_dXY[comboIdx][cIdx]->Fill(heli[0]); 
                        hists.fake_unused_pt2_MD0_dZ[comboIdx][cIdx]->Fill(heli[1]);
                        hists.fake_unused_pt2_MD1_dXY[comboIdx][cIdx]->Fill(heli[2]); 
                        hists.fake_unused_pt2_MD1_dZ[comboIdx][cIdx]->Fill(heli[3]);
                    }

                    if (rz_simple.first > -900) {
                        hists.fake_unused_pt2_MD0_rz_simple[comboIdx][cIdx]->Fill(rz_simple.first);
                        hists.fake_unused_pt2_MD1_rz_simple[comboIdx][cIdx]->Fill(rz_simple.second);
                    }
                }
            }      
        } 
    }
         if (writeRoot) {
            // Check for duplicates
            for (size_t i = 0; i < pT2_matched_simIdx.size(); i++) {
                bool isDup = false;
                if (!pT2_matched_simIdx[i].empty()) {
                    int simIdx = pT2_matched_simIdx[i][0];
                    // If this sim track was found more than once, it's a duplicate
                    if (simIdx >= 0 && simIdx < sim_pT2_matched.size() && sim_pT2_matched[simIdx] > 1) {
                        isDup = true;
                    }
                }
                pT2_isDuplicate.push_back(isDup ? 1 : 0);
            }
                        // === INJECT pT2 INTO TC ===
            if (!reader.tc_pt || !reader.sim_tcIdx || !reader.sim_tcIdxAll) {
                std::cerr << "ERROR: tc_ or sim_tc branches are not loaded in rootReader!" << std::endl;
                return 1; 
            }

            for (size_t i = 0; i < pT2_pt.size(); i++) {
                
                // 1. Get the new global ID this track will have in the main array
                int new_tc_index = reader.tc_pt->size();

                // 2. Inject the track
                reader.tc_pt->push_back(pT2_pt[i]);
                reader.tc_eta->push_back(pT2_eta[i]);
                reader.tc_phi->push_back(pT2_phi[i]);
                reader.tc_isFake->push_back(pT2_isFake[i]);
                reader.tc_isDuplicate->push_back(pT2_isDuplicate[i]);
                reader.tc_type->push_back(2); // Spoof as pLS to pass the filter
                
                int primary_simIdx = pT2_matched_simIdx[i].empty() ? -1 : pT2_matched_simIdx[i][0];
                reader.tc_simIdx->push_back(primary_simIdx);
                reader.tc_simIdxAll->push_back(pT2_matched_simIdx[i]);
                reader.tc_nhits->push_back(2);
                reader.tc_nlayers->push_back(2);
                
                // pT2 has 0 Outer Tracker hits
                if (reader.tc_nhitOT) reader.tc_nhitOT->push_back(0);

                // If fake, match frac is 0. If real, let's say 1.0 (100% matched)
                if (reader.tc_pMatched) reader.tc_pMatched->push_back(pT2_isFake[i] ? 0.0 : 1.0);

                // tc_simIdxAllFrac must be the same length as the tc_simIdxAll array
                if (reader.tc_simIdxAllFrac) {
                    std::vector<float> frac;
                    if (!pT2_matched_simIdx[i].empty()) {
                        frac.push_back(1.0); // 100% match
                    }
                    reader.tc_simIdxAllFrac->push_back(frac);
                }

                // Indices to higher level LST objects. Since pT2 doesn't use these, they get -1 (Null)
                if (reader.tc_pt5Idx) reader.tc_pt5Idx->push_back(-1);
                if (reader.tc_pt3Idx) reader.tc_pt3Idx->push_back(-1);
                if (reader.tc_t5Idx)  reader.tc_t5Idx->push_back(-1);

                // We DO know the underlying pLS index, so we can save it!
                if (reader.tc_plsIdx) reader.tc_plsIdx->push_back(pT2_plsIdx[i]);

                // 3. === THE REVERSE TRUTH MATCH ===
                // Tell the simulated particle that this new track rescued it!
                if (primary_simIdx >= 0 && primary_simIdx < reader.sim_tcIdxAll->size()) {
                    // Add our new track ID to the particle's list of successful matches
                    reader.sim_tcIdxAll->at(primary_simIdx).push_back(new_tc_index);
                    
                    // If the particle was completely lost before, officially mark it as found
                    if (reader.sim_tcIdx->at(primary_simIdx) == -1) {
                        reader.sim_tcIdx->at(primary_simIdx) = new_tc_index;
                        //std::cout << "found" <<std::endl;
                    }
                   // reader.sim_tcIdx->at(primary_simIdx) = new_tc_index;
                }
            }
            // ==========================
            // CRITICAL: Fill the TTree entry for this event!
            outTree->Fill();
            trainTree->Fill();
        }
    }
   // --- AUTOMATED CUT CALCULATION ---
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "IDEAL CUTS PER CATEGORY (Target Efficiency: " << targetPercent << "%)" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    for (int i = 0; i < 13; ++i) {
        for (int c = 0; c < 2; ++c) { // Loop over charge (0 = Pos, 1 = Neg)

            // Skip if this specific geometry+charge combo has no events
            if (hists.real_pt2_deltaPT[i][c]->GetEntries() == 0) continue;

            // 1-Sided configuration (upper bound only)
            double q_1[1];
            double p_1[1] = { targetPercent / 100.0 };

            // 2-Sided configuration (symmetric tails)
            double tail = (1.0 - (targetPercent / 100.0)) / 2.0;
            double q_2[2];
            double p_2[2] = { tail, 1.0 - tail };

            // --- Extrapolation Variables (1-Sided) ---
            double cut_dZ0 = 0, cut_dZ1 = 0;
            double cut_dXY0 = 0, cut_dXY1 = 0;

            // --- Simple R-Z & Kinematic Variables (2-Sided) ---
            double cut_RZ0Min = 0, cut_RZ0Max = 0;
            double cut_RZ1Min = 0, cut_RZ1Max = 0;
            double cut_dPhiMin = 0, cut_dPhiMax = 0;
            double cut_dPtMin = 0,  cut_dPtMax = 0;

            // --- LST Variables (2-Sided) ---
            double cut_LSTdPhiMin = 0,    cut_LSTdPhiMax = 0;
            double cut_LSTdBetaMin = 0,   cut_LSTdBetaMax = 0;
            double cut_LSTbetaOutMin = 0, cut_LSTbetaOutMax = 0;
            double cut_LSTOrgZMin = 0,    cut_LSTOrgZMax = 0;
            double cut_LSTKinZMin = 0,    cut_LSTKinZMax = 0;

            // --- Calculate 1-Sided Quantiles ---
            if (hists.real_pt2_MD0_dZ[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD0_dZ[i][c]->GetQuantiles(1, q_1, p_1);
                cut_dZ0 = q_1[0];
            }
            if (hists.real_pt2_MD1_dZ[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD1_dZ[i][c]->GetQuantiles(1, q_1, p_1);
                cut_dZ1 = q_1[0];
            }
            if (hists.real_pt2_MD0_dXY[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD0_dXY[i][c]->GetQuantiles(1, q_1, p_1);
                cut_dXY0 = q_1[0];
            }
            if (hists.real_pt2_MD1_dXY[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD1_dXY[i][c]->GetQuantiles(1, q_1, p_1);
                cut_dXY1 = q_1[0];
            }

            // --- Calculate 2-Sided Quantiles (Standard) ---
            if (hists.real_pt2_MD0_rz_simple[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD0_rz_simple[i][c]->GetQuantiles(2, q_2, p_2);
                cut_RZ0Min = q_2[0]; cut_RZ0Max = q_2[1];
            }
            if (hists.real_pt2_MD1_rz_simple[i][c]->GetEntries() > 0) {
                hists.real_pt2_MD1_rz_simple[i][c]->GetQuantiles(2, q_2, p_2);
                cut_RZ1Min = q_2[0]; cut_RZ1Max = q_2[1];
            }
            if (hists.real_pt2_deltaPHI[i][c]->GetEntries() > 0) {
                hists.real_pt2_deltaPHI[i][c]->GetQuantiles(2, q_2, p_2);
                cut_dPhiMin = q_2[0]; cut_dPhiMax = q_2[1];
            }
            if (hists.real_pt2_deltaPT[i][c]->GetEntries() > 0) {
                hists.real_pt2_deltaPT[i][c]->GetQuantiles(2, q_2, p_2);
                cut_dPtMin = q_2[0]; cut_dPtMax = q_2[1];
            }

            // --- Calculate 2-Sided Quantiles (LST) ---
            if (hists.real_pt2_LSTdPhi[i][c]->GetEntries() > 0) {
                hists.real_pt2_LSTdPhi[i][c]->GetQuantiles(2, q_2, p_2);
                cut_LSTdPhiMin = q_2[0]; cut_LSTdPhiMax = q_2[1];
            }
            if (hists.real_pt2_LSTdBeta[i][c]->GetEntries() > 0) {
                hists.real_pt2_LSTdBeta[i][c]->GetQuantiles(2, q_2, p_2);
                cut_LSTdBetaMin = q_2[0]; cut_LSTdBetaMax = q_2[1];
            }
            if (hists.real_pt2_LSTbetaOut[i][c]->GetEntries() > 0) {
                hists.real_pt2_LSTbetaOut[i][c]->GetQuantiles(2, q_2, p_2);
                cut_LSTbetaOutMin = q_2[0]; cut_LSTbetaOutMax = q_2[1];
            }
            if (hists.real_pt2_LSTOrgZRes[i][c]->GetEntries() > 0) {
                hists.real_pt2_LSTOrgZRes[i][c]->GetQuantiles(2, q_2, p_2);
                cut_LSTOrgZMin = q_2[0]; cut_LSTOrgZMax = q_2[1];
            }
            if (hists.real_pt2_LSTKinZRes[i][c]->GetEntries() > 0) {
                hists.real_pt2_LSTKinZRes[i][c]->GetQuantiles(2, q_2, p_2);
                cut_LSTKinZMin = q_2[0]; cut_LSTKinZMax = q_2[1];
            }

            // --- Print Formatting ---
            std::cout << "\n>>> CATEGORY: " << hists.catTitles[i] << " | CHARGE: " << hists.chargeTitles[c] << " <<<" << std::endl;
            std::cout << std::fixed << std::setprecision(4);

            std::cout << "  [Helical Extrapolation Cuts]" << std::endl;
            std::cout << "    MD0 dZ Cut:          < " << cut_dZ0 << " cm" << std::endl;
            std::cout << "    MD1 dZ Cut:          < " << cut_dZ1 << " cm" << std::endl;
            std::cout << "    MD0 dXY Cut:         < " << cut_dXY0 << " cm" << std::endl;
            std::cout << "    MD1 dXY Cut:         < " << cut_dXY1 << " cm" << std::endl;

            std::cout << "[Simple Pointing Cuts]" << std::endl;
            std::cout << "    MD0 R-Z Res:           " << std::setw(8) << cut_RZ0Min << " to " << cut_RZ0Max << " cm" << std::endl;
            std::cout << "    MD1 R-Z Res:           " << std::setw(8) << cut_RZ1Min << " to " << cut_RZ1Max << " cm" << std::endl;

            std::cout << "  [Kinematic Cuts]" << std::endl;
            std::cout << "    Delta Phi:             " << std::setw(8) << cut_dPhiMin << " to " << cut_dPhiMax << " rad" << std::endl;
            std::cout << "    Delta pT:              " << std::setw(8) << cut_dPtMin << " to " << cut_dPtMax << " GeV" << std::endl;

            std::cout << "  [LST Component Cuts]" << std::endl;
            std::cout << "    LST Delta Phi:         " << std::setw(8) << cut_LSTdPhiMin << " to " << cut_LSTdPhiMax << " rad" << std::endl;
            std::cout << "    LST Delta Beta:        " << std::setw(8) << cut_LSTdBetaMin << " to " << cut_LSTdBetaMax << " rad" << std::endl;
            std::cout << "    LST Beta Out:          " << std::setw(8) << cut_LSTbetaOutMin << " to " << cut_LSTbetaOutMax << " rad" << std::endl;
            std::cout << "    LST Geometric Z-Res:   " << std::setw(8) << cut_LSTOrgZMin << " to " << cut_LSTOrgZMax << " cm" << std::endl;
            std::cout << "    LST Kinematic Z-Res:   " << std::setw(8) << cut_LSTKinZMin << " to " << cut_LSTKinZMax << " cm" << std::endl;
        }
    }
    std::cout << "\n" << std::string(80, '=') << "\n" << std::endl;

    if (makePlots){
        auto recipes = getPt2Recipes(hists);
        Plotting plotter; 
        plotter.plotRecipes(recipes, outputDir);
    }
    
    // === NEW: Write to Output Root File ===
    if (writeRoot) {
        outRootFile->cd();
        outTree->Write(); // Write the TTree containing pT3 branches + new pT2 branches!
        outRootFile->Close();
        std::cout << "Saved ROOT file containing old + new pT2 branches to: " << outputDir << "/LSTNtuple_with_pT2.root" << std::endl;
        trainFile->cd();
        trainTree->Write();
        trainFile->Close();
        std::cout << "Saved training data to: " << trainFileName << "\n";
        }
    // =======================================


    return 0;
}