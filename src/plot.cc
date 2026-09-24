#include "plot.h"

#include "histograms.h"
#include "plot_recipes.h"
#include "plotting.h"

int runPlot(const Config& cfg)
{
    HistogramManager hists;
    hists.load(cfg.histFile);

    Plotting plotter;
    plotter.plotRecipes(getPt2Recipes(hists), cfg.outputDir);
    return 0;
}
