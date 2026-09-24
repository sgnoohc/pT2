#ifndef CONFIG_H
#define CONFIG_H

#include <string>

enum class Mode { Process, Scan, Plot };

struct Config {
    Mode mode = Mode::Process;

    // process
    bool lowPT = false;
    bool writeRoot = false;
    int nEvents = -1;
    std::string inputFile;       // defaults depend on lowPT
    std::string pixelMapDir;     // defaults depend on lowPT
    std::string nnModelDir;      // holds model.onnx, mean.npy, std.npy

    // scan
    double targetPercent = 90;

    // all modes
    std::string outputDir = "output";
    std::string histFile;        // defaults to <outputDir>/pt2_hists.root
};

// Parse "pt2 <process|scan|plot> [options]". Returns false (after printing usage) on bad input.
bool parseArgs(int argc, char** argv, Config& cfg);

void printConfig(const Config& cfg);

#endif
