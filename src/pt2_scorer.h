#ifndef PT2_SCORER_H
#define PT2_SCORER_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <onnxruntime_cxx_api.h>

#include "pt2.h"
#include "rootReader.h"

// ONNX NN that scores a pT2 candidate as real (-> 1) or fake (-> 0)
class Pt2Scorer
{
public:
    static constexpr size_t kNFeatures = 21;
    using FeatureVector = std::array<float, kNFeatures>;

    // meanPath / stdPath: .npy normalization arrays saved by the training script
    Pt2Scorer(const std::string &modelPath, const std::string &meanPath, const std::string &stdPath);

    // Model inputs, in the order of FEATURE_NAMES in pt2_ml/train_*.py
    static FeatureVector features(const rootReader &reader, const pT2 &pt2);

    // Score many candidates in as few ONNX calls as possible; one score per row, same order
    std::vector<float> score(const std::vector<FeatureVector> &raw) const;
    float score(const FeatureVector &raw) const;
    float score(const rootReader &reader, const pT2 &pt2) const;

private:
    Ort::Env env_;
    Ort::SessionOptions sessionOptions_;
    std::unique_ptr<Ort::Session> session_;
    Ort::MemoryInfo memInfo_;

    std::vector<float> mean_, std_;
    std::string inputName_, outputName_;
    bool dynamicBatch_ = false; // model accepts N rows per call (else 1 row per call)

    // Run the model on n already-normalized rows
    void run(float *normed, size_t n, float *out) const;
};

#endif
