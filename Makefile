#Makefile for INCA Fortran project

.DEFAULT_GOAL := inca.exe

SRC_DIR = src
OBJ_DIR = obj

FC = gfortran

# --- CONFIGURATION: CHOOSE DEBUG OR RELEASE ---
# 1. DEBUG MODE (Slow, use for fixing errors)
#FFLAGS = -fopenmp -Wall -g -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -J$(OBJ_DIR)

# 2. RELEASE MODE (Fast, use for actual calculation)
# We use -O3 for maximum speed and -fopenmp for parallelization
FFLAGS = -fopenmp -O3 -march=native -funroll-loops -J$(OBJ_DIR)

# Library linking (LAPACK/BLAS)
LIB = -llapack -lblas

SRCF90 = $(wildcard $(SRC_DIR)/*.f90)
SRCF   = $(wildcard $(SRC_DIR)/*.f)
OBJF90 = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCF90))
OBJF   = $(patsubst $(SRC_DIR)/%.f,$(OBJ_DIR)/%.o,$(SRCF))
OBJ    = $(OBJF90) $(OBJF)

$(shell mkdir -p $(OBJ_DIR))
# Pattern rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_DIR) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -I$(OBJ_DIR) -c $< -o $@
# Final executable
inca.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB)
# Generate dependencies
deps:
	python3 generate_deps.py
# Clean
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod inca.exe test.log
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod 
# Test
test:
	python3 tests/run_tests.py 
# Dependencies (generated)
include deps.mk
