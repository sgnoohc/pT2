#include "tools.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include "rootReader.h"
#include <vector>
#include <algorithm>
#include <array>

// Some Useful Constants
constexpr int superbin_nEta = 25;
constexpr int superbin_nPhi = 72;
constexpr int superbin_nZ   = 25;
constexpr float superbin_etaMin = -2.6f;
constexpr float superbin_etaMax =  2.6f;
constexpr float superbin_zMin   = -30.0f;
constexpr float superbin_zMax   =  30.0f;
constexpr float superbin_ptBins_norm[2] = {0.8f, 2.0f};
constexpr float superbin_ptBins_low[2] = {0.6f, 2.0f};


float deltaPhi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    if (dphi >  M_PI)
        dphi -= 2.0f * M_PI;
    else if (dphi < -M_PI)
        dphi += 2.0f * M_PI;
    return dphi;
}


float deltaEta(float eta1, float eta2) {
    return eta1 - eta2;
}


float deltaPt(float pt1, float pt2) {
    return pt1 - pt2;
}


std::vector<int> getDetIdsForLS(const rootReader& reader, size_t k) {
    std::vector<int> detIds;
    if (!reader.ls_mdIdx0 || !reader.ls_mdIdx1 || !reader.md_detId)
        throw std::runtime_error("Required LS or MD branches not found");
    if (k >= reader.ls_mdIdx0->size())
        throw std::out_of_range("LS index out of range");
    int mdIdx0 = reader.ls_mdIdx0->at(k);
    if (mdIdx0 >= 0 && mdIdx0 < (int)reader.md_detId->size()) {
        detIds.push_back(reader.md_detId->at(mdIdx0));
    }
    return detIds;
}

static void loadSingleSuperbinDetIdMap(const std::string& filename, SuperbinToDetIdMap& superbinToDetIds)
{
    superbinToDetIds.clear();
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open pixelMap file: " + filename);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        int superbin = -1;
        int nDet = 0;
        if (!(ss >> superbin >> nDet)) {
            throw std::runtime_error("Failed to read pixelMap at line: " + line);
        }
        if (nDet == 0)
            continue;
        std::vector<int> detIds;
        detIds.reserve(nDet);
        for (int i = 0; i < nDet; ++i) {
            int detId;
            if (!(ss >> detId)) {
                throw std::runtime_error("Expected more detIDs in pixelMap at line: " + line);
            }
            detIds.push_back(detId);
        }
        superbinToDetIds[superbin] = std::move(detIds);
    }
}

void loadSuperbinDetIdMap(const std::string& dir, SuperbinToDetIdMap& pos, SuperbinToDetIdMap& neg, SuperbinToDetIdMap& non)
{
    std::string posFile = dir + "pLS_map_pos_ElCheapo.txt";
    std::string negFile = dir + "pLS_map_neg_ElCheapo.txt";
    std::string nonFile = dir + "pLS_map_ElCheapo.txt";

    loadSingleSuperbinDetIdMap(posFile, pos);
    loadSingleSuperbinDetIdMap(negFile, neg);
    loadSingleSuperbinDetIdMap(nonFile, non);
}

int CalculateSuperbin(const rootReader& reader, size_t j, bool lowPT) {
    float pt = reader.pls_pt->at(j);
    float phi = reader.pls_phi->at(j);
    float eta = reader.pls_eta->at(j);
    float z = reader.pls_origin_z.at(j);
    // pT bin
    const auto& ptBins = (lowPT ? superbin_ptBins_low : superbin_ptBins_norm);
    int ipt = (pt < ptBins[1]) ? 0 : 1;
    // eta bin 
    int ieta;
    if (eta < superbin_etaMin) {
        ieta = 0;
    } else if (eta > superbin_etaMax) {
        ieta = superbin_nEta - 1;
    } else {
        const float etaWidth = (superbin_etaMax - superbin_etaMin) / superbin_nEta;
        ieta = static_cast<int>((eta - superbin_etaMin) / etaWidth);
        if (ieta < 0) ieta = 0;
        if (ieta >= superbin_nEta) ieta = superbin_nEta - 1;
    }
    // phi bin
    while (phi < -M_PI) phi += 2.0f * M_PI;
    while (phi >  M_PI) phi -= 2.0f * M_PI;
    const float phiWidth = (2.0f * M_PI) / superbin_nPhi;
    int iphi = static_cast<int>((phi + M_PI) / phiWidth);
    if (iphi < 0) iphi = 0;
    if (iphi >= superbin_nPhi) iphi = superbin_nPhi - 1;
    // z bin
    int iz;
    if (z < superbin_zMin) {
        iz = 0;
    } else if (z > superbin_zMax) {
        iz = superbin_nZ - 1;
    } else {
        const float zWidth = (superbin_zMax - superbin_zMin) / superbin_nZ;
        iz = static_cast<int>((z - superbin_zMin) / zWidth);
        if (iz < 0) iz = 0;
        if (iz >= superbin_nZ) iz = superbin_nZ - 1;
    }
    // compute superbin 
    const int superbin =
        (ipt  * superbin_nPhi * superbin_nEta * superbin_nZ) +
        (ieta * superbin_nPhi * superbin_nZ) +
        (iphi * superbin_nZ) +
        iz;
    return superbin;
}


/*float CalculatePlsZ(const rootReader& reader, size_t j) {
    //TODO: This is just a rough approximation that I don't really understand
    int nhit = reader.pls_nhit->at(j);
    if (nhit != 3 && nhit != 4) {
        throw std::runtime_error("CalculatePlsZ: nhit must be 3 or 4");
    }
    // Grab the endpoints for now 
    std::array<float,3> p0 = { reader.pls_hit0_x->at(j), reader.pls_hit0_y->at(j), reader.pls_hit0_z->at(j) };
    std::array<float,3> p1;
    if (nhit == 3) {
        p1 = { reader.pls_hit2_x->at(j), reader.pls_hit2_y->at(j), reader.pls_hit2_z->at(j) };
    } else { 
        p1 = { reader.pls_hit3_x->at(j), reader.pls_hit3_y->at(j), reader.pls_hit3_z->at(j) };
    }
    float dx = p1[0] - p0[0];
    float dy = p1[1] - p0[1];
    float dz = p1[2] - p0[2];
    // --- Solve for t where x=y=0 ---
    float tx = (dx != 0.0f) ? -p0[0] / dx : 0.0f;
    float ty = (dy != 0.0f) ? -p0[1] / dy : 0.0f;
    float t = 0.0f;
    int count = 0;
    if (dx != 0.0f) { t += tx; ++count; }
    if (dy != 0.0f) { t += ty; ++count; }
    if (count == 0) return p0[2]; 
    t /= count;
    return p0[2] + t * dz;
}*/

float CalculatePlsZ(const rootReader& reader, size_t j) {
    int nhit = reader.pls_nhit->at(j);
    if (nhit != 3 && nhit != 4) {
        throw std::runtime_error("CalculatePlsZ: nhit must be 3 or 4");
    }

    // 1. Get the coordinates of the inner-most hit (hit 0)
    float x0 = reader.pls_hit0_x->at(j);
    float y0 = reader.pls_hit0_y->at(j);
    float z0 = reader.pls_hit0_z->at(j);

    // 2. Get the coordinates of the outer-most hit (hit 2 or 3)
    float x1 = (nhit == 3) ? reader.pls_hit2_x->at(j) : reader.pls_hit3_x->at(j);
    float y1 = (nhit == 3) ? reader.pls_hit2_y->at(j) : reader.pls_hit3_y->at(j);
    float z1 = (nhit == 3) ? reader.pls_hit2_z->at(j) : reader.pls_hit3_z->at(j);

    // 3. Calculate the transverse radius R from the beamline
    float r0 = std::sqrt(x0*x0 + y0*y0);
    float r1 = std::sqrt(x1*x1 + y1*y1);

    // 4. Calculate the changes in R and Z
    float dr = r1 - r0;
    float dz = z1 - z0;

    // Protect against division by zero (e.g., if a track somehow has two hits at the exact same radius)
    if (dr == 0.0f) {
        return z0;
    }

    // 5. Linear extrapolation in the R-Z plane back to the beamline (R = 0)
    // Formula: z_origin = z0 - r0 * (slope)
    return z0 - r0 * (dz / dr);
}


void buildPt2sForPLS(size_t pls_idx, 
                    const rootReader& reader, 
                    const SuperbinToDetIdMap& superbinToDetIds_POS, 
                    const SuperbinToDetIdMap& superbinToDetIds_NEG, 
                    const SuperbinToDetIdMap& superbinToDetIds_NON, 
                    const DetIdToLSMap& detIdToLS, 
                    pT2Collection& pt2s)
{
    // =========================================================
    // SANITY CHECK MODE: No pixel map / superbin cuts applied.
    // Pair this pLS with EVERY single LS in the event.
    // =========================================================

    size_t nLS = reader.ls_pt->size(); // Get the total number of LS
    
    for (size_t ls_idx = 0; ls_idx < nLS; ++ls_idx) {
        if (pt2TruthFinder(reader, pls_idx, ls_idx)) {
        pT2 obj(pls_idx, ls_idx);
        pt2s.push_back(obj); }
    }
/*    int superbin = reader.pls_superbin[pls_idx];
    // Determine which SuperbinToDetIdMap to use
    const SuperbinToDetIdMap* selectedMap = nullptr;
    float pt = reader.pls_pt->at(pls_idx);
    int charge = reader.pls_charge->at(pls_idx);
    if (pt >= 2.0f) {
        selectedMap = &superbinToDetIds_NON;
    } 
    else if (charge > 0) {
        selectedMap = &superbinToDetIds_POS;
    } 
    else if (charge < 0) {
        selectedMap = &superbinToDetIds_NEG;
    } 
    else {
        throw std::runtime_error("pLS has charge = 0");
    }
    // Lookup detIds in the selected map
    auto itDetIds = selectedMap->find(superbin);
    if (itDetIds == selectedMap->end()) return;
    const std::vector<int>& detIds = itDetIds->second;
    for (int detId : detIds) {
        auto itLS = detIdToLS.find(detId);
        if (itLS == detIdToLS.end()) continue;
        const std::vector<size_t>& lsIndices = itLS->second;
        for (size_t ls_idx : lsIndices) {
            pT2 obj(pls_idx, ls_idx);
            pt2s.push_back(obj);
        }
    }*/
}


bool pt2TruthFinder(const rootReader& reader, size_t plsIdx, size_t lsIdx)
{
    int plsSimIdx = reader.pls_simIdx->at(plsIdx);
    int lsSimIdx  = reader.ls_simIdx->at(lsIdx);
    if (plsSimIdx != lsSimIdx) return false;
    const auto& plsFracVec = reader.pls_simIdxAllFrac->at(plsIdx);
    float plsFrac  = plsFracVec.empty()  ? 0.0f : plsFracVec[0];
    int plsNhit   = reader.pls_nhit->at(plsIdx);
    const auto& lsFracVec = reader.ls_simIdxAllFrac->at(lsIdx);
    float lsFrac  = lsFracVec.empty()  ? 0.0f : lsFracVec[0];
    int lsNhit    = 4;  
    float realFraction = ((plsFrac * plsNhit) + (lsFrac * lsNhit)) / (plsNhit + lsNhit);
    return (realFraction > 0.75f);
}


UsedMask buildUsedMask(const rootReader& reader) {
    size_t nLS  = reader.ls_pt  ? reader.ls_pt->size()  : 0;
    size_t nPLS = reader.pls_pt ? reader.pls_pt->size() : 0;
    std::vector<bool> ls_isUsed(nLS, false);
    std::vector<bool> pls_isUsed(nPLS, false);
    // Helper function to mark an index as used
    auto markIndex = [](std::vector<bool>& mask, int idx) {
        if (idx >= 0 && idx < (int)mask.size())
            mask[idx] = true;
    };
    // TC_PT5
    if (reader.tc_pt5Idx && reader.pt5_plsIdx && reader.pt5_t5Idx) {
        for (int tcIdx : *reader.tc_pt5Idx) {
            if (tcIdx < 0 || tcIdx >= (int)reader.pt5_plsIdx->size()) continue;
            int plsIdx = reader.pt5_plsIdx->at(tcIdx);
            int t5Idx  = reader.pt5_t5Idx->at(tcIdx);
            markIndex(pls_isUsed, plsIdx);
            // T5 -> T3 -> LS
            if (t5Idx >= 0) {
                if (reader.t5_t3Idx0 && t5Idx < (int)reader.t5_t3Idx0->size()) {
                    int t3Idx0 = reader.t5_t3Idx0->at(t5Idx);
                    if (t3Idx0 >= 0) {
                        if (reader.t3_lsIdx0) markIndex(ls_isUsed, reader.t3_lsIdx0->at(t3Idx0));
                        if (reader.t3_lsIdx1) markIndex(ls_isUsed, reader.t3_lsIdx1->at(t3Idx0));
                    }
                }
                if (reader.t5_t3Idx1 && t5Idx < (int)reader.t5_t3Idx1->size()) {
                    int t3Idx1 = reader.t5_t3Idx1->at(t5Idx);
                    if (t3Idx1 >= 0) {
                        if (reader.t3_lsIdx0) markIndex(ls_isUsed, reader.t3_lsIdx0->at(t3Idx1));
                        if (reader.t3_lsIdx1) markIndex(ls_isUsed, reader.t3_lsIdx1->at(t3Idx1));
                    }
                }
            }
        }
    }
    // TC_T5
    if (reader.tc_t5Idx) {
        for (int t5Idx : *reader.tc_t5Idx) {
            if (t5Idx < 0) continue;
            if (reader.t5_t3Idx0 && t5Idx < (int)reader.t5_t3Idx0->size()) {
                int t3Idx0 = reader.t5_t3Idx0->at(t5Idx);
                if (t3Idx0 >= 0) {
                    if (reader.t3_lsIdx0) markIndex(ls_isUsed, reader.t3_lsIdx0->at(t3Idx0));
                    if (reader.t3_lsIdx1) markIndex(ls_isUsed, reader.t3_lsIdx1->at(t3Idx0));
                }
            }
            if (reader.t5_t3Idx1 && t5Idx < (int)reader.t5_t3Idx1->size()) {
                int t3Idx1 = reader.t5_t3Idx1->at(t5Idx);
                if (t3Idx1 >= 0) {
                    if (reader.t3_lsIdx0) markIndex(ls_isUsed, reader.t3_lsIdx0->at(t3Idx1));
                    if (reader.t3_lsIdx1) markIndex(ls_isUsed, reader.t3_lsIdx1->at(t3Idx1));
                }
            }
        }
    }
    // TC_PT3
    if (reader.tc_pt3Idx && reader.pt3_plsIdx && reader.pt3_t3Idx) {
        for (int pt3Idx : *reader.tc_pt3Idx) {
            if (pt3Idx < 0 || pt3Idx >= (int)reader.pt3_plsIdx->size()) continue;
            int plsIdx = reader.pt3_plsIdx->at(pt3Idx);
            markIndex(pls_isUsed, plsIdx);
            int t3Idx = reader.pt3_t3Idx->at(pt3Idx);
            if (t3Idx >= 0) {
                if (reader.t3_lsIdx0) markIndex(ls_isUsed, reader.t3_lsIdx0->at(t3Idx));
                if (reader.t3_lsIdx1) markIndex(ls_isUsed, reader.t3_lsIdx1->at(t3Idx));
            }
        }
    }
    return {ls_isUsed, pls_isUsed};
}


bool pt2UsedCalculator(const rootReader& reader, size_t plsIdx, size_t lsIdx) {
    bool plsUsed = reader.pls_isUsed[plsIdx];
    bool lsUsed  = reader.ls_isUsed[lsIdx];
    return plsUsed || lsUsed;
}
