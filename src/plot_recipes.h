#ifndef PLOT_RECIPES_H
#define PLOT_RECIPES_H

#include <string>
#include <vector>
#include "histograms.h"
#include <TH1.h>

struct PlotRecipe {
    std::string title;
    std::string xAxis;
    std::string yAxis;
    std::string filename;
    std::vector<TH1D*> hists; 
    std::vector<std::string> legend{};
    bool logX = false;
    bool logY = false;
    bool normalize = false;
    bool printYields = false;
};

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists);

#endif
