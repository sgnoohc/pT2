#ifndef TOOLS_H
#define TOOLS_H

#include <unordered_map>
#include <vector>
#include <string>
#include "rootReader.h"
#include "pt2.h"

// Type aliases
using DetIdToLSMap = std::unordered_map<int, std::vector<size_t>>;
using SuperbinToDetIdMap = std::unordered_map<int, std::vector<int>>;

float deltaPhi(float phi1, float phi2);
float deltaEta(float eta1, float eta2);
float deltaPt(float pt1, float pt2);

float CalculatePlsZ(const rootReader& reader, size_t j);
int CalculateSuperbin(const rootReader& reader, size_t j, bool lowPT);

std::vector<int> getDetIdsForLS(const rootReader& reader, size_t k);

struct UsedMask {
    std::vector<bool> ls_isUsed;
    std::vector<bool> pls_isUsed;
};

UsedMask buildUsedMask(const rootReader& reader);

// superbin -> detID maps (POS / NEG / NON)
void loadSuperbinDetIdMap(
    const std::string& dir,
    SuperbinToDetIdMap& pos,
    SuperbinToDetIdMap& neg,
    SuperbinToDetIdMap& non
);

bool pt2TruthFinder(const rootReader& reader, size_t plsIdx, size_t lsIdx);

bool pt2UsedCalculator(const rootReader& reader, size_t plsIdx, size_t lsIdx);

// pt2 builder
void buildPt2sForPLS(size_t pls_idx,
                     const rootReader& reader,
                     const SuperbinToDetIdMap& superbinToDetIds_POS,
                     const SuperbinToDetIdMap& superbinToDetIds_NEG,
                     const SuperbinToDetIdMap& superbinToDetIds_NON,
                     const DetIdToLSMap& detIdToLS,
                     pT2Collection& pt2s);

#endif
