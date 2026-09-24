#ifndef PT2_NTUPLE_WRITER_H
#define PT2_NTUPLE_WRITER_H

#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "rootReader.h"
#include "pt2.h"

// Clones the input LST ntuple, adds pT2 branches, and injects pT2s into the tc_ collections.
class Pt2NtupleWriter {
public:
    Pt2NtupleWriter(const std::string& path, rootReader& reader);
    ~Pt2NtupleWriter();

    // Call after reader.GetEntry()
    void beginEvent();
    // Record one selected pT2
    void add(const pT2& pt2);
    // Flag duplicates, inject pT2s as TCs, and fill the tree
    void endEvent();
    // Write the tree and close the file
    void close();

    const std::string& path() const { return path_; }

private:
    void flagDuplicates();
    void injectIntoTCs();

    std::string path_;
    rootReader& reader_;
    TFile* file_ = nullptr;
    TTree* tree_ = nullptr;

    std::vector<float> pt_, eta_, phi_;
    std::vector<int> plsIdx_, lsIdx_;
    std::vector<int> isFake_, isUsed_, isDuplicate_;
    std::vector<float> deltaPt_, deltaEta_, deltaPhi_, dR_;
    std::vector<float> nnScore_;
    std::vector<std::vector<int>> matchedSimIdx_;
    std::vector<int> simMatched_;  // per sim track: number of pT2s matched to it
};

#endif
