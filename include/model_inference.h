#pragma once

#include <onnxruntime_cxx_api.h>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <array>

// =============================================================================
// Tiny .npy loader — supports float32 and float64, 1D arrays only
// =============================================================================
inline std::vector<float> loadNpy(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open: " + path);

    char magic[6];
    f.read(magic, 6);
    if (std::strncmp(magic, "\x93NUMPY", 6) != 0)
        throw std::runtime_error("Not a .npy file: " + path);

    uint8_t major, minor;
    f.read(reinterpret_cast<char*>(&major), 1);
    f.read(reinterpret_cast<char*>(&minor), 1);

    uint32_t hlen = 0;
    if (major == 1) {
        uint16_t h16;
        f.read(reinterpret_cast<char*>(&h16), 2);
        hlen = h16;
    } else {
        f.read(reinterpret_cast<char*>(&hlen), 4);
    }

    std::string header(hlen, ' ');
    f.read(&header[0], hlen);

    bool isDouble = (header.find("'f8'") != std::string::npos ||
                     header.find("float64") != std::string::npos ||
                     header.find("<f8") != std::string::npos);

    size_t shapePos = header.find("shape");
    size_t lp = header.find('(', shapePos);
    size_t rp = header.find(')', lp);

    std::string shapeStr = header.substr(lp + 1, rp - lp - 1);
    while (!shapeStr.empty() && (shapeStr.back() == ',' || shapeStr.back() == ' '))
        shapeStr.pop_back();

    size_t n = std::stoul(shapeStr);

    std::vector<float> out(n);

    if (isDouble) {
        std::vector<double> tmp(n);
        f.read(reinterpret_cast<char*>(tmp.data()), n * sizeof(double));
        for (size_t i = 0; i < n; ++i)
            out[i] = static_cast<float>(tmp[i]);
    } else {
        f.read(reinterpret_cast<char*>(out.data()), n * sizeof(float));
    }

    return out;
}

// =============================================================================
// Pt2Scorer — ONNX inference wrapper
// =============================================================================
class Pt2Scorer {
public:
    static constexpr size_t N_FEATURES = 21;

    Pt2Scorer(const std::string& modelPath,
              const std::string& meanPath,
              const std::string& stdPath)
        : env_(ORT_LOGGING_LEVEL_WARNING, "pt2scorer"),
          sessionOptions_(),
          memInfo_(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault))
    {
        sessionOptions_.SetIntraOpNumThreads(1);
        sessionOptions_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

        session_ = std::make_unique<Ort::Session>(env_, modelPath.c_str(), sessionOptions_);

        // Load normalization
        mean_ = loadNpy(meanPath);
        std_  = loadNpy(stdPath);

        if (mean_.size() != N_FEATURES || std_.size() != N_FEATURES)
            throw std::runtime_error("mean/std size mismatch with N_FEATURES");

        // Get input/output names SAFELY
        Ort::AllocatorWithDefaultOptions allocator;

        auto inputNameAlloc  = session_->GetInputNameAllocated(0, allocator);
        auto outputNameAlloc = session_->GetOutputNameAllocated(0, allocator);

        inputName_  = inputNameAlloc.release();
        outputName_ = outputNameAlloc.release();

        // Optional: print once for sanity
        std::cout << "ONNX input: " << inputName_
                  << " | output: " << outputName_ << std::endl;
    }

    ~Pt2Scorer() {
        Ort::AllocatorWithDefaultOptions allocator;
        if (inputName_)  allocator.Free(inputName_);
        if (outputName_) allocator.Free(outputName_);
    }

    // Score a single candidate
    float score(const std::array<float, N_FEATURES>& raw) const {
        std::array<float, N_FEATURES> normed;

        // Normalize + protect against NaN/inf
        for (size_t i = 0; i < N_FEATURES; ++i) {
            float val = raw[i];
            if (!std::isfinite(val)) val = 0.0f;
            normed[i] = (val - mean_[i]) / std_[i];
        }

        // Create tensor
        const int64_t shape[2] = {1, static_cast<int64_t>(N_FEATURES)};

        Ort::Value inputTensor = Ort::Value::CreateTensor<float>(
            memInfo_,
            const_cast<float*>(normed.data()),
            N_FEATURES,
            shape,
            2
        );

        const char* inputNames[]  = {inputName_};
        const char* outputNames[] = {outputName_};

        auto outputs = session_->Run(
            Ort::RunOptions{nullptr},
            inputNames,  &inputTensor, 1,
            outputNames, 1
        );

        // Safety check
        if (outputs.empty() || !outputs[0].IsTensor())
            throw std::runtime_error("Invalid ONNX output");

        float logit = outputs[0].GetTensorMutableData<float>()[0];

        // Sigmoid
        return 1.0f / (1.0f + std::exp(-logit));
    }

private:
    Ort::Env env_;
    Ort::SessionOptions sessionOptions_;
    std::unique_ptr<Ort::Session> session_;

    Ort::MemoryInfo memInfo_;

    std::vector<float> mean_;
    std::vector<float> std_;

    char* inputName_  = nullptr;
    char* outputName_ = nullptr;
};
