#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>

#include "rootReader.h"
#include "extrapolation.h"

using namespace std;

int main(int argc, char** argv) {
    string inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple.root";
    string outputFile = "nr_convergence_study.root";
    int maxEvents = 1000; // Keep it small for a quick study

    cout << "Starting Newton-Raphson Convergence Study..." << endl;

    rootReader reader;
    if (!reader.Init(inputFile, "tree")) {
        cerr << "Failed to open input file!" << endl;
        return 1;
    }

    // ==========================================
    // INITIALIZE PLOTS
    // ==========================================
    TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");

    // TProfiles show the *average* value on the Y-axis for each X-bin (iteration)
    TProfile* prof_dist_vs_iter = new TProfile("prof_dist_vs_iter", 
        "Average 3D Distance vs Iteration;Iteration Step;Avg Distance to Target [cm]", 15, 0, 15);
    
    TProfile* prof_step_vs_iter = new TProfile("prof_step_vs_iter", 
        "Average Step Size (#Delta#phi) vs Iteration;Iteration Step;Avg |#Delta#phi| [rad]", 15, 0, 15);

    TH1D* h_iterations_needed = new TH1D("h_iterations_needed", 
        "Iterations Required to Converge (#Delta#phi < 1e-6);Number of Iterations;Track Count", 15, 0, 15);

    Long64_t entries = reader.GetEntries();
    if (maxEvents > 0 && maxEvents < entries) entries = maxEvents;

    // ==========================================
    // MAIN EVENT LOOP
    // ==========================================
    for (Long64_t ievt = 0; ievt < entries; ++ievt) {
        reader.GetEntry(ievt);
        if (ievt % 100 == 0) cout << "Processing event " << ievt << " / " << entries << endl;

        size_t nPLS = reader.pls_pt->size();
        size_t nLS = reader.ls_pt->size();

        // Loop over pLS and LS to find TRUE matching tracks (Signal)
        for (size_t p = 0; p < nPLS; ++p) {
            int pls_sim = reader.pls_simIdx->at(p);
            if (pls_sim < 0) continue; // Ignore fakes/noise

            for (size_t l = 0; l < nLS; ++l) {
                int ls_sim = reader.ls_simIdx->at(l);
                
                // Only run study on REAL physical track pairs
                if (pls_sim == ls_sim) {
                    
                    // 1. EXTRACT DATA FOR HELIX (Matching extrapolation.cpp logic)
                    vector<ROOT::Math::XYVector> hits_xy;
                    vector<double> hits_z;
                    
                    if (reader.pls_hit0_x->at(p) > -900) { hits_xy.emplace_back(reader.pls_hit0_x->at(p), reader.pls_hit0_y->at(p)); hits_z.push_back(reader.pls_hit0_z->at(p)); }
                    if (reader.pls_hit1_x->at(p) > -900) { hits_xy.emplace_back(reader.pls_hit1_x->at(p), reader.pls_hit1_y->at(p)); hits_z.push_back(reader.pls_hit1_z->at(p)); }
                    if (reader.pls_hit2_x->at(p) > -900) { hits_xy.emplace_back(reader.pls_hit2_x->at(p), reader.pls_hit2_y->at(p)); hits_z.push_back(reader.pls_hit2_z->at(p)); }
                    if (reader.pls_hit3_x->at(p) > -900) { hits_xy.emplace_back(reader.pls_hit3_x->at(p), reader.pls_hit3_y->at(p)); hits_z.push_back(reader.pls_hit3_z->at(p)); }

                    if (hits_xy.size() < 2) continue;

                    double R_pred_cm = (reader.pls_pt->at(p) / (0.3 * 3.8)) * 100.0;
                    ROOT::Math::XYVector center = extrapolation::fitCircleWithFixedRadius(hits_xy, R_pred_cm);
                    if (center.X() == 0 && center.Y() == 0) continue;

                    // Pitch logic
                    double first_angle = TMath::ATan2(hits_xy.front().Y() - center.Y(), hits_xy.front().X() - center.X());
                    double last_angle = TMath::ATan2(hits_xy.back().Y() - center.Y(), hits_xy.back().X() - center.X());
                    while (last_angle - first_angle > TMath::Pi()) last_angle -= 2 * TMath::Pi();
                    while (last_angle - first_angle < -TMath::Pi()) last_angle += 2 * TMath::Pi();
                    
                    double dPhi = last_angle - first_angle;
                    if (std::abs(dPhi) < 1e-9) continue;

                    double a = (hits_z.back() - hits_z.front()) / dPhi;
                    double b = hits_z.front() - a * first_angle;

                    // 2. GET TARGET MD (We will just use MD0 for the study)
                    int md0 = reader.ls_mdIdx0->at(l);
                    if (md0 < 0) continue;
                    double tx = reader.md_anchor_x->at(md0);
                    double ty = reader.md_anchor_y->at(md0);
                    double tz = reader.md_anchor_z->at(md0);

                    // 3. THE NEWTON-RAPHSON STUDY LOOP
                    double phi0 = std::atan2(ty - center.Y(), tx - center.X());
                    while (phi0 - last_angle >  TMath::Pi()) phi0 -= 2*TMath::Pi();
                    while (phi0 - last_angle < -TMath::Pi()) phi0 += 2*TMath::Pi();

                    double p_fine = phi0;
                    bool converged = false;
                    int iters_taken = 15;

                    for (int iter = 0; iter < 15; ++iter) {
                        double c = std::cos(p_fine), s = std::sin(p_fine);
                        
                        // Calculate current 3D physical distance
                        double hX = center.X() + R_pred_cm * c;
                        double hY = center.Y() + R_pred_cm * s;
                        double hZ = a * p_fine + b;
                        double current_dist = std::sqrt(std::pow(tx - hX, 2) + std::pow(ty - hY, 2) + std::pow(tz - hZ, 2));

                        // Record Distance for this iteration
                        prof_dist_vs_iter->Fill(iter, current_dist);

                        // Calculate Newton Step
                        double f = (tx - center.X() - R_pred_cm*c)*(R_pred_cm*s) + 
                                   (ty - center.Y() - R_pred_cm*s)*(-R_pred_cm*c) + 
                                   a*(a*p_fine + b - tz);
                        double fp = R_pred_cm*R_pred_cm + a*a;
                        
                        double step = f / fp;
                        
                        // Record Step Size for this iteration
                        prof_step_vs_iter->Fill(iter, std::abs(step));

                        p_fine -= step;

                        // Check Convergence (Tolerance = 1 micro-radian)
                        if (std::abs(step) < 1e-6 && !converged) {
                            iters_taken = iter + 1;
                            converged = true;
                        }
                    }
                    
                    h_iterations_needed->Fill(iters_taken);
                }
            }
        }
    }

// ==========================================
    // SAVE RESULTS & DRAW PLOTS (Updated Styling)
    // ==========================================
    outFile->cd();
    prof_dist_vs_iter->Write();
    prof_step_vs_iter->Write();
    h_iterations_needed->Write();

    // Create PNG visuals
    gStyle->SetOptStat(0);
    TCanvas c1("c1", "Convergence", 1200, 400);
    c1.Divide(3, 1);

     // --- Plot 1: Distance vs Iteration ---
    c1.cd(1);
    TGraph* g_dist = new TGraph();
    int pt_idx1 = 0;
    for(int i = 1; i <= prof_dist_vs_iter->GetNbinsX(); ++i) {
        // Only plot points that actually have data
        if(prof_dist_vs_iter->GetBinEntries(i) > 0) {
            g_dist->SetPoint(pt_idx1++, prof_dist_vs_iter->GetBinCenter(i), prof_dist_vs_iter->GetBinContent(i));
        }
    }
    g_dist->SetTitle("Average 3D Distance vs Iteration;Iteration Step;Avg Distance to Target [cm]");
    g_dist->SetMarkerStyle(20);          
    g_dist->SetMarkerColor(kBlue+2);     
    g_dist->SetMarkerSize(1.2);
    g_dist->SetLineColor(kBlue-4);       
    g_dist->SetLineStyle(2);             
    g_dist->SetLineWidth(2);
    // "A" = draw Axes, "P" = draw Points, "L" = draw Line connecting them
    g_dist->Draw("APL");                  

    // --- Plot 2: Step Size vs Iteration ---
    c1.cd(2);
    TGraph* g_step = new TGraph();
    int pt_idx2 = 0;
    for(int i = 1; i <= prof_step_vs_iter->GetNbinsX(); ++i) {
        if(prof_step_vs_iter->GetBinEntries(i) > 0) {
            g_step->SetPoint(pt_idx2++, prof_step_vs_iter->GetBinCenter(i), prof_step_vs_iter->GetBinContent(i));
        }
    }
    g_step->SetTitle("Average Step Size (#Delta#phi) vs Iteration;Iteration Step;Avg |#Delta#phi| [rad]");
    g_step->SetMarkerStyle(20);          
    g_step->SetMarkerColor(kRed+1);      
    g_step->SetMarkerSize(1.2);
    g_step->SetLineColor(kRed-4);        
    g_step->SetLineStyle(2);             
    g_step->SetLineWidth(2);
    g_step->Draw("APL");
    gPad->SetLogy();                           // Keep Log scale to see the exponential drop

    // --- Plot 3: Iterations Needed ---
    c1.cd(3);
    h_iterations_needed->SetFillColor(kGreen-5);    // Softer green fill
    h_iterations_needed->SetLineColor(kGreen+2);    // Darker green outline
    h_iterations_needed->SetLineWidth(2);
    h_iterations_needed->Draw("HIST");

    // Save exactly as a PNG!
    c1.SaveAs("NewtonRaphson_Convergence.png");

    outFile->Close();
    cout << "Done! Results saved to nr_convergence_study.root and NewtonRaphson_Convergence.png" << endl;
    return 0;
}
