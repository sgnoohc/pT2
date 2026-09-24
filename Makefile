CXX      := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra
ONNX_DIR := /cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ENV_FILES/onnxruntime-linux-x64-1.17.0

CPPFLAGS := -MMD -MP $(shell root-config --cflags) -isystem $(ONNX_DIR)/include
LDLIBS   := $(shell root-config --libs) -L$(ONNX_DIR)/lib -lonnxruntime -Wl,-rpath,$(ONNX_DIR)/lib

SRC_DIR   := src
BUILD_DIR := build
BIN_DIR   := bin

# Files with their own main(); everything else in src/ is shared library code
MAINS     := main study_convergence
LIB_SRCS  := $(filter-out $(addprefix $(SRC_DIR)/,$(addsuffix .cc,$(MAINS))),$(wildcard $(SRC_DIR)/*.cc))
LIB_OBJS  := $(LIB_SRCS:$(SRC_DIR)/%.cc=$(BUILD_DIR)/%.o)

all: $(BIN_DIR)/pt2

$(BIN_DIR)/pt2: $(BUILD_DIR)/main.o $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Optional: make study_convergence
study_convergence: $(BIN_DIR)/study_convergence
$(BIN_DIR)/study_convergence: $(BUILD_DIR)/study_convergence.o $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(wildcard $(BUILD_DIR)/*.d)

.PHONY: all clean study_convergence
