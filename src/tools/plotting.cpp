#include "plotting.h"
#include "plot_recipes.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TH1.h>
#include <THStack.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <fstream>

namespace fs = std::filesystem;

//-------------------------------------------------------------
// Core plotting loop
//-------------------------------------------------------------
void Plotting::plotRecipes(const std::vector<PlotRecipe>& recipes,
                           const std::string& outputDir)
{
    if (!fs::exists(outputDir))
        fs::create_directories(outputDir);

    gStyle->SetOptStat(0);

    int colors[] = {kSpring+10, kRed-7, kGreen-6, kMagenta-7, kRed-4, kCyan-6};

    for (const auto& recipe : recipes)
    {
        // Check histograms exist
        if (std::any_of(recipe.hists.begin(), recipe.hists.end(),
                        [](TH1* h){ return !h; }))
        {
            std::cerr << "Warning: Missing histogram for plot \""
                      << recipe.title << "\"\n";
            continue;
        }

        bool useLegend = !recipe.legend.empty() &&
                         recipe.legend.size() == recipe.hists.size();

        for (bool useLogY : {false, recipe.logY})
        {
            //--------------------------------------------------
            // Canvas
            //--------------------------------------------------
            TCanvas canvas(("canvas_" + recipe.filename).c_str(),
                           recipe.title.c_str(), 900, 700);

            //--------------------------------------------------
            // Pads
            //--------------------------------------------------
            TPad* mainPad;
            TPad* yieldPad = nullptr;

            if (recipe.printYields)
            {
                // MAIN PAD (top 75%)
                mainPad = new TPad("mainPad","",0,0.25,1,1);
                mainPad->SetTopMargin(0.08);
                mainPad->SetBottomMargin(0.15);
                mainPad->SetLeftMargin(0.12);
                mainPad->SetRightMargin(0.05);
                mainPad->Draw();

                // YIELD PAD (bottom 25%)
                yieldPad = new TPad("yieldPad","",0,0,1,0.25);
                yieldPad->SetTopMargin(0.08);
                yieldPad->SetBottomMargin(0.28);
                yieldPad->SetLeftMargin(0.12);
                yieldPad->SetRightMargin(0.05);
                yieldPad->Draw();
            }
            else
            {
                mainPad = &canvas;
                mainPad->SetTopMargin(0.08);
                mainPad->SetBottomMargin(0.17);
                mainPad->SetLeftMargin(0.12);
                mainPad->SetRightMargin(0.05);
            }

            //--------------------------------------------------
            // Prepare histograms
            //--------------------------------------------------
            mainPad->cd();
            mainPad->SetLogy(useLogY);

            std::vector<TH1D*> hists;
            double ymax = 0.0;

            for (size_t i = 0; i < recipe.hists.size(); ++i)
            {
                TH1D* h = (TH1D*)recipe.hists[i]->Clone(
                    (std::string(recipe.hists[i]->GetName()) + "_" + recipe.filename).c_str()
                );
                h->SetDirectory(nullptr);

                if (recipe.normalize)
                {
                    double integral = h->Integral();
                    if (integral > 0) h->Scale(1.0 / integral);

                    h->SetFillColorAlpha(colors[i % 6], 0.35);
                    h->SetLineColor(colors[i % 6]);
                    h->SetLineWidth(3);
                }
                else
                {
                    h->SetFillColor(colors[i % 6]);
                    h->SetLineColor(kBlack);
                    h->SetLineWidth(1);
                    h->SetFillStyle(1001);
                }

                ymax = std::max(ymax, h->GetMaximum());
                hists.push_back(h);
            }

            //--------------------------------------------------
            // Draw histograms
            //--------------------------------------------------
            if (recipe.normalize)
            {
                for (size_t i = 0; i < hists.size(); ++i)
                {
                    TH1D* h = hists[i];
                    h->SetMaximum(useLogY ? ymax*20.0 : ymax*1.35);

                    if (i == 0) h->Draw("HIST");
                    else h->Draw("HIST SAME");

                    h->GetXaxis()->SetTitle(recipe.xAxis.c_str());
                    h->GetXaxis()->SetTitleSize(0.06);
                    h->GetXaxis()->SetLabelSize(0.05);
                    h->GetXaxis()->SetTitleOffset(1.0);

                    h->GetYaxis()->SetTitle(recipe.yAxis.c_str());
                    h->GetYaxis()->SetTitleSize(0.06);
                    h->GetYaxis()->SetLabelSize(0.05);
                    h->GetYaxis()->SetTitleOffset(0.9);
                }
            }
            else
            {
                THStack* stack = new THStack(("stack_" + recipe.filename).c_str(), recipe.title.c_str());
                for (auto* h : hists) stack->Add(h);

                stack->Draw("HIST");
                double stackMax = stack->GetMaximum();
                stack->SetMaximum(useLogY ? stackMax*5.0 : stackMax*1.25);
                if (useLogY) stack->SetMinimum(0.001);

                stack->GetXaxis()->SetTitle(recipe.xAxis.c_str());
                stack->GetXaxis()->SetTitleSize(0.06);
                stack->GetXaxis()->SetLabelSize(0.05);
                stack->GetXaxis()->SetTitleOffset(1.1);

                stack->GetYaxis()->SetTitle(recipe.yAxis.c_str());
                stack->GetYaxis()->SetTitleSize(0.06);
                stack->GetYaxis()->SetLabelSize(0.05);
                stack->GetYaxis()->SetTitleOffset(1.0);
            }

            //--------------------------------------------------
            // Legend
            //--------------------------------------------------
            if (useLegend)
            {
                TLegend* legend = new TLegend(0.68, 0.72, 0.90, 0.90);
               // legend->SetBorderSize(0);
               // legend->SetFillStyle(0);
                legend->SetTextSize(0.045);

                for (size_t i = 0; i < hists.size(); ++i)
                    legend->AddEntry(hists[i], recipe.legend[i].c_str(),
                                     recipe.normalize ? "l" : "f");

                legend->Draw();
            }

            //--------------------------------------------------
            // YIELD PAD (now counts underflow + overflow)
            //--------------------------------------------------
            if (recipe.printYields && yieldPad)
            {
                yieldPad->cd();
                yieldPad->Clear();

                double x = 0.05;
                double y = 0.80;
                double yStep = 0.16;

                for (size_t i = 0; i < hists.size(); ++i)
                {
                    TH1D* h = hists[i];

                    // ---- FIX: include underflow + overflow ----
                    int nbins = h->GetNbinsX();
                    double yield = h->Integral(0, nbins + 1);

                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextColor(h->GetLineColor());
                    latex.SetTextSize(0.16);

                    latex.DrawLatex(x, y,
                        (recipe.legend[i] + " : " +
                        std::to_string(static_cast<long long>(yield))).c_str());

                    y -= yStep;
                }
            }

            //--------------------------------------------------
            // Save
            //--------------------------------------------------
            canvas.cd();
            canvas.Modified();
            canvas.Update();

            std::string suffix = useLogY ? "_logY" : "";
            canvas.SaveAs((outputDir + "/" + recipe.filename + suffix + ".png").c_str());
            canvas.SaveAs((outputDir + "/" + recipe.filename + suffix + ".pdf").c_str());
        }
    }

    std::cout << "All plots saved to " << outputDir << "\n";
    generateHTMLGallery(outputDir);
}
//-------------------------------------------------------------
// HTML gallery generator
//-------------------------------------------------------------
void Plotting::generateHTMLGallery(const std::string& outputDir)
{
    std::vector<fs::path> pngFiles;
    for (const auto& entry : fs::directory_iterator(outputDir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".png") {
            pngFiles.push_back(entry.path().filename());
        }
    }

    std::sort(pngFiles.begin(), pngFiles.end());

    std::ofstream html(outputDir + "/index.html");
    if (!html.is_open()) {
        std::cerr << "Error: Cannot create HTML file in " << outputDir << "\n";
        return;
    }

    html << R"(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Plot Gallery</title>
<style>
body { font-family: Arial, sans-serif; margin: 20px; }
input { padding: 5px; margin-bottom: 10px; width: 300px; }
#zoom { margin-left: 10px; }
.plot-container { display: flex; flex-wrap: wrap; gap: 20px; }
.plot-item { text-align: center; }
.plot-item img { width: 200px; transition: width 0.2s; cursor: pointer; }
</style>
</head>
<body>

<h1>Plot Gallery</h1>

<input type="text" id="search" placeholder="Search plots...">
<label>Zoom: <input type="range" id="zoom" min="100" max="500" value="200"></label>

<div class="plot-container">
)";

    for (const auto& pngPath : pngFiles) {
        std::string name = pngPath.stem().string();

        html << "<div class=\"plot-item\" data-title=\"" << name << "\">"
             << "<a href=\"" << name << ".pdf\" target=\"_blank\">"
             << "<img src=\"" << pngPath.string() << "\" alt=\"" << name << "\">"
             << "</a>"
             << "<div>" << name << "</div>"
             << "</div>\n";
    }

    html << R"(</div>
<script>
window.onload = function() {
    const searchInput = document.getElementById('search');
    const zoomSlider = document.getElementById('zoom');
    const items = document.querySelectorAll('.plot-item');

    searchInput.addEventListener('input', () => {
        const query = searchInput.value.toLowerCase();
        items.forEach(item => {
            const title = item.getAttribute('data-title').toLowerCase();
            item.style.display = title.includes(query) ? 'block' : 'none';
        });
    });

    zoomSlider.addEventListener('input', () => {
        const size = zoomSlider.value + 'px';
        items.forEach(item => {
            item.querySelector('img').style.width = size;
        });
    });
};
</script>
</body>
</html>)";

    html.close();
    std::cout << "HTML gallery generated at " << outputDir << "/index.html\n";
}
