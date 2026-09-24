#ifndef PROCESS_H
#define PROCESS_H

#include "config.h"

// Event loop: build pT2s, fill histograms into cfg.histFile, optionally write the pT2 ntuple
int runProcess(const Config& cfg);

#endif
