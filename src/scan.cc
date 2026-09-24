#include "scan.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

#include "histograms.h"

namespace {

// Upper bound keeping `eff` of the entries (0 if empty)
double upperCut(TH1D* h, double eff)
{
    if (h->GetEntries() == 0) return 0;
    double q[1], p[1] = {eff};
    h->GetQuantiles(1, q, p);
    return q[0];
}

// Symmetric-tail window keeping `eff` of the entries ({0, 0} if empty)
std::pair<double, double> windowCut(TH1D* h, double eff)
{
    if (h->GetEntries() == 0) return {0, 0};
    double tail = (1.0 - eff) / 2.0;
    double q[2], p[2] = {tail, 1.0 - tail};
    h->GetQuantiles(2, q, p);
    return {q[0], q[1]};
}

void printUpper(const std::string& label, double cut, const std::string& unit)
{
    std::cout << "    " << std::left << std::setw(22) << label << std::right << "< " << cut << " " << unit << "\n";
}

void printWindow(const std::string& label, std::pair<double, double> cut, const std::string& unit)
{
    std::cout << "    " << std::left << std::setw(22) << label << std::right
              << std::setw(8) << cut.first << " to " << cut.second << " " << unit << "\n";
}

} // namespace

int runScan(const Config& cfg)
{
    HistogramManager hists;
    hists.load(cfg.histFile);

    double eff = cfg.targetPercent / 100.0;

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "IDEAL CUTS PER CATEGORY (Target Efficiency: " << cfg.targetPercent << "%)\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << std::fixed << std::setprecision(4);

    for (int i = 0; i < kNCat; ++i) {
        for (int c = 0; c < kNCharge; ++c) {
            const Pt2HistSet& h = hists.set(true, false, i, c);
            if (h.deltaPT->GetEntries() == 0) continue;

            std::cout << "\n>>> CATEGORY: " << hists.catTitles[i] << " | CHARGE: " << hists.chargeTitles[c] << " <<<\n";

            std::cout << "  [Helical Extrapolation Cuts]\n";
            printUpper("MD0 dZ Cut:", upperCut(h.MD0_dZ, eff), "cm");
            printUpper("MD1 dZ Cut:", upperCut(h.MD1_dZ, eff), "cm");
            printUpper("MD0 dXY Cut:", upperCut(h.MD0_dXY, eff), "cm");
            printUpper("MD1 dXY Cut:", upperCut(h.MD1_dXY, eff), "cm");

            std::cout << "  [Simple Pointing Cuts]\n";
            printWindow("MD0 R-Z Res:", windowCut(h.MD0_rz_simple, eff), "cm");
            printWindow("MD1 R-Z Res:", windowCut(h.MD1_rz_simple, eff), "cm");

            std::cout << "  [Kinematic Cuts]\n";
            printWindow("Delta Phi:", windowCut(h.deltaPHI, eff), "rad");
            printWindow("Delta pT:", windowCut(h.deltaPT, eff), "GeV");

            std::cout << "  [LST Component Cuts]\n";
            printWindow("LST Delta Phi:", windowCut(h.LSTdPhi, eff), "rad");
            printWindow("LST Delta Beta:", windowCut(h.LSTdBeta, eff), "rad");
            printWindow("LST Beta Out:", windowCut(h.LSTbetaOut, eff), "rad");
            printWindow("LST Geometric Z-Res:", windowCut(h.LSTOrgZRes, eff), "cm");
            printWindow("LST Kinematic Z-Res:", windowCut(h.LSTKinZRes, eff), "cm");
        }
    }
    std::cout << "\n" << std::string(80, '=') << "\n" << std::endl;

    return 0;
}
