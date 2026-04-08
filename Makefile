# SPLASH Fortran MPI Library Makefile
# This Makefile builds the SPLASH static library only.
# It does NOT build examples by default.

.NOTPARALLEL:

# Compiler settings
FC = mpif90
FFLAGS = -O3 -g

# Directory structure
SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include
LIB_DIR = lib

# Library name
LIB_NAME = libsplash.a

# Installation prefix (configurable by user)
# Safe default: install inside current project directory
# Override on command line for system-style install, e.g.:
# make install PREFIX=/usr/local
PREFIX ?= $(CURDIR)/install

# SPLASH source files in exact dependency order (DO NOT reorder)
SPLASH_SRCS = \
  $(SRC_DIR)/SPLASH_Precision.f90 \
  $(SRC_DIR)/SPLASH_Parameters.f90 \
  $(SRC_DIR)/SPLASH_MPI_Constants.f90 \
  $(SRC_DIR)/SPLASH_Buffer.f90 \
  $(SRC_DIR)/SPLASH_FFT_Pre.f90 \
  $(SRC_DIR)/SPLASH_LT_Tuning.f90 \
  $(SRC_DIR)/SPLASH_MPI_Part.f90 \
  $(SRC_DIR)/SPLASH_FFT_kernel.f90 \
  $(SRC_DIR)/SPLASH_LocalTranspose.f90 \
  $(SRC_DIR)/SPLASH_GlobalTranspose.f90 \
  $(SRC_DIR)/SPLASH_Poisson_Algebra.f90 \
  $(SRC_DIR)/SPLASH_3D_FFT.f90 \
  $(SRC_DIR)/SPLASH_Repartitioning.f90 \
  $(SRC_DIR)/SPLASH_Halos_Exchange.f90 \
  $(SRC_DIR)/SPLASH_3D_Periodic_Poisson.f90 \
  $(SRC_DIR)/SPLASH_Integrated_Solver.f90 \
  $(SRC_DIR)/SPLASH_Planner.f90 \
  $(SRC_DIR)/SPLASH.f90

# Object files in build directory, preserving source order
SPLASH_OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SPLASH_SRCS))

.PHONY: all install clean print-prefix

# Default target: build the library only
all: $(LIB_DIR)/$(LIB_NAME)

# Optional helper: show current install prefix
print-prefix:
	@echo "PREFIX = $(PREFIX)"

# Create directories if needed
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(INCLUDE_DIR):
	mkdir -p $(INCLUDE_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# Compile source files to object files in build/, modules in include/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR) $(INCLUDE_DIR)
	$(FC) $(FFLAGS) -I$(INCLUDE_DIR) -J$(INCLUDE_DIR) -c $< -o $@

# Archive object files into static library
$(LIB_DIR)/$(LIB_NAME): $(SPLASH_OBJS) | $(LIB_DIR)
	ar rcs $@ $^

# Install library and module files to PREFIX
install: all
	mkdir -p $(PREFIX)/lib $(PREFIX)/include
	cp $(LIB_DIR)/$(LIB_NAME) $(PREFIX)/lib/
	cp $(INCLUDE_DIR)/*.mod $(PREFIX)/include/
	@echo "Installed library to: $(PREFIX)/lib/$(LIB_NAME)"
	@echo "Installed module files to: $(PREFIX)/include/"

# Clean intermediate and built artifacts
clean:
	rm -rf $(BUILD_DIR) $(INCLUDE_DIR) $(LIB_DIR) 