# -----------------------------------------------------------------
# Makefile for BayesOmics
# Supported platforms:
#   Unix / Linux (UNIX)
#   MacOS (MAC)
# -----------------------------------------------------------------

OUTPUT = BayesOmics64-static

# Check the operating system using uname
SYS := $(shell uname -s)

ifeq ($(SYS),Linux)
	# Path to dependent libraries
	EIGEN       = $(EIGEN3_INCLUDE_DIR)
	BOOST_ROOT  = $(BOOST_ROOT_DIR)
	BOOST_LIB   = $(BOOST_LIB_DIR)
	ARMA        = $(ARMA_INCLUDE_DIR)
	HDF5_INCLUDE = $(HDF5_INCLUDE_DIR)
	HDF5_LIB     = $(HDF5_LIB_DIR)
endif

ifeq ($(SYS),MAC)
	EIGEN       = ~/.local/lib/eigen-3.4.0
	BOOST_ROOT  = ~/.local/lib/boost_1_80_0
	ARMA        = /opt/homebrew/Cellar/armadillo/12.6.2/include
	CC          = /opt/homebrew/opt/llvm/bin/clang
	CXX_UNIX    = /opt/homebrew/opt/llvm/bin/clang++
	LDFLAGS     = -L/opt/homebrew/opt/llvm/lib -Wl,-rpath,/opt/homebrew/opt/llvm/lib -L/opt/homebrew/lib
	CXXFLAGS    = "-I/opt/homebrew/opt/llvm/include -I/opt/homebrew/opt/llvm/include/c++/v1 -I/opt/homebrew/include"
	AR          = /opt/homebrew/opt/llvm/bin/llvm-ar
	LD          = /opt/homebrew/opt/llvm/bin/clang
endif

# -----------------------------------------------------------------
# System specific flags and settings
ifeq ($(SYS),Linux)
	CXX_LINUX = g++
	CXXFLAGS  = -fopenmp -w -O3 -m64 -DEIGEN_NO_DEBUG -msse2 -std=c++17 -I.
	CXXFLAGS += -I $(EIGEN) -I $(BOOST_ROOT) -I $(ARMA) -I$(HDF5_INCLUDE) -DUNIX -m64 -DBOOST_NO_CXX11_SCOPED_ENUMS
	LDFLAGS  += -L$(BOOST_LIB) -L$(HDF5_LIB)

	# Link HDF5 first
	LDLIBS   += -lhdf5 -lhdf5_hl -lsz -laec

	# Boost (depends on zlib)
	LDLIBS   += -lboost_iostreams -lboost_system -lboost_filesystem -lboost_program_options -lboost_regex

	# Compression (after Boost)
	LDLIBS   += -lz

	# Math + threading
	LDLIBS   += -static -DARMA_DONT_USE_WRAPPER -lopenblas -lgfortran -lquadmath -lgomp -lpthread -ldl

	CXX       = $(CXX_LINUX)
endif

ifeq ($(SYS),MAC)
	CXX_LINUX = g++
	CXXFLAGS  = -g -fopenmp -w -O3 -m64 -DEIGEN_NO_DEBUG -msse2 -std=c++17 -I.
	CXXFLAGS += -I $(EIGEN) -I $(BOOST_ROOT) -I $(ARMA) -DUNIX -m64 -ldl
	LDFLAGS  += -L$(BOOST_LIB_DIR)/ -lhdf5_cpp -lhdf5
	LDLIBS   += -lboost_filesystem -lboost_iostreams -lboost_system -lboost_program_options -lboost_regex
	LDLIBS   += /opt/homebrew/Cellar/armadillo/12.6.2/lib/libarmadillo.dylib
	CXX       = $(CXX_UNIX)
endif

# -----------------------------------------------------------------
# Directory settings
SRC_DIR     = ../src
BIN_DIR     = ./bin
INCLUDE_DIR = ../include

CXXFLAGS   += -I $(INCLUDE_DIR)

# Source and object files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)
TARGET  = $(BIN_DIR)/${OUTPUT}

# -----------------------------------------------------------------
# Git commit hash
GIT_COMMIT_HASH := $(shell git rev-parse HEAD)

# -----------------------------------------------------------------
# Main target
all: clean $(TARGET)

# Linking the target executable
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compiling source files with the Git commit hash
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -DGIT_COMMIT_HASH=\"$(GIT_COMMIT_HASH)\" -c -o $@ $<

# Creating necessary directories
$(BIN_DIR):
	mkdir -p $@

# Cleaning up the build
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR)

# Phony targets
.PHONY: all clean
