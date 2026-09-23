CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -MMD -MP

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

SRC_DIR   := src
TOOLS_DIR := $(SRC_DIR)/tools
INC_DIR   := include
BUILD_DIR := build
BIN_DIR   := bin
BIN       := $(BIN_DIR)/pt2

SOURCES := $(SRC_DIR)/main.cpp $(wildcard $(TOOLS_DIR)/*.cpp)
OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(SOURCES)))

all: $(BIN)

$(BIN): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(ROOTLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(TOOLS_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(OBJECTS:.o=.d)

.PHONY: all clean
