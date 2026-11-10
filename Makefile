FC = mpif90
AR = ar
FFLAGS = -O3
# FFLAGS += -ffree-line-length-none -fallow-argument-mismatch

SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include
LIB_DIR = lib
LIB_NAME = libisfft.a
PREFIX ?= $(CURDIR)/install

.NOTPARALLEL:
.DELETE_ON_ERROR:

ISFFT_SRCS = \
  $(SRC_DIR)/ISFFT_Communicator.f90 \
  $(SRC_DIR)/ISFFT_Precision.f90 \
  $(SRC_DIR)/ISFFT_Problem_State.f90 \
  $(SRC_DIR)/ISFFT_Discrete_Operators.f90 \
  $(SRC_DIR)/ISFFT_Workspace.f90 \
  $(SRC_DIR)/ISFFT_Plans.f90 \
  $(SRC_DIR)/ISFFT_MPI_Part.f90 \
  $(SRC_DIR)/ISFFT_Halos_Exchange.f90 \
  $(SRC_DIR)/ISFFT_FFT_kernel.f90 \
  $(SRC_DIR)/ISFFT_DST_kernel.f90 \
  $(SRC_DIR)/ISFFT_DCT_kernel.f90 \
  $(SRC_DIR)/ISFFT_Modal_Operations.f90 \
  $(SRC_DIR)/ISFFT_Communication.f90 \
  $(SRC_DIR)/ISFFT_Axis_Execution.f90 \
  $(SRC_DIR)/ISFFT_Transform_Execution.f90 \
  $(SRC_DIR)/ISFFT_Allocation.f90 \
  $(SRC_DIR)/ISFFT_Precompute.f90 \
  $(SRC_DIR)/ISFFT_Planning.f90 \
  $(SRC_DIR)/ISFFT_Initialization.f90 \
  $(SRC_DIR)/ISFFT_Redistribution.f90 \
  $(SRC_DIR)/ISFFT_Integrated_Solver.f90 \
  $(SRC_DIR)/ISFFT.f90

ISFFT_OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(ISFFT_SRCS))

.PHONY: all library install clean print-prefix
all: library
library: $(LIB_DIR)/$(LIB_NAME)
	@printf 'ISFFT library is ready.\n'

# Recover when generated module files have been removed.
ifneq ($(words $(wildcard $(INCLUDE_DIR)/isfft*.mod)),$(words $(ISFFT_OBJS)))
.PHONY: missing-modules
$(ISFFT_OBJS): missing-modules
endif

$(BUILD_DIR) $(INCLUDE_DIR) $(LIB_DIR):
	@mkdir -p $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 Makefile | $(BUILD_DIR) $(INCLUDE_DIR)
	@$(FC) $(FFLAGS) -I$(INCLUDE_DIR) -J$(INCLUDE_DIR) -c $< -o $@

# Each object depends on the modules imported by its source file.
$(BUILD_DIR)/ISFFT_Problem_State.o: $(BUILD_DIR)/ISFFT_Precision.o
$(BUILD_DIR)/ISFFT_Discrete_Operators.o: $(BUILD_DIR)/ISFFT_Precision.o
$(BUILD_DIR)/ISFFT_Workspace.o: $(BUILD_DIR)/ISFFT_Precision.o
$(BUILD_DIR)/ISFFT_Plans.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Problem_State.o
$(BUILD_DIR)/ISFFT_MPI_Part.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Problem_State.o
$(BUILD_DIR)/ISFFT_Halos_Exchange.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o $(BUILD_DIR)/ISFFT_MPI_Part.o
$(BUILD_DIR)/ISFFT_FFT_kernel.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Workspace.o
$(BUILD_DIR)/ISFFT_DST_kernel.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_FFT_kernel.o
$(BUILD_DIR)/ISFFT_DCT_kernel.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_FFT_kernel.o
$(BUILD_DIR)/ISFFT_Modal_Operations.o: $(BUILD_DIR)/ISFFT_Workspace.o
$(BUILD_DIR)/ISFFT_Communication.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o
$(BUILD_DIR)/ISFFT_Axis_Execution.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Workspace.o $(BUILD_DIR)/ISFFT_MPI_Part.o \
  $(BUILD_DIR)/ISFFT_DCT_kernel.o $(BUILD_DIR)/ISFFT_DST_kernel.o $(BUILD_DIR)/ISFFT_FFT_kernel.o
$(BUILD_DIR)/ISFFT_Transform_Execution.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_Axis_Execution.o $(BUILD_DIR)/ISFFT_Communication.o
$(BUILD_DIR)/ISFFT_Allocation.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_Communication.o $(BUILD_DIR)/ISFFT_Transform_Execution.o
$(BUILD_DIR)/ISFFT_Precompute.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Discrete_Operators.o $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o \
  $(BUILD_DIR)/ISFFT_Workspace.o $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_FFT_kernel.o \
  $(BUILD_DIR)/ISFFT_Allocation.o
$(BUILD_DIR)/ISFFT_Planning.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_DCT_kernel.o $(BUILD_DIR)/ISFFT_DST_kernel.o \
  $(BUILD_DIR)/ISFFT_FFT_kernel.o $(BUILD_DIR)/ISFFT_Axis_Execution.o $(BUILD_DIR)/ISFFT_Communication.o \
  $(BUILD_DIR)/ISFFT_Transform_Execution.o $(BUILD_DIR)/ISFFT_Allocation.o $(BUILD_DIR)/ISFFT_Precompute.o
$(BUILD_DIR)/ISFFT_Initialization.o: $(BUILD_DIR)/ISFFT_Communicator.o $(BUILD_DIR)/ISFFT_Precision.o \
  $(BUILD_DIR)/ISFFT_Plans.o $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_Communication.o $(BUILD_DIR)/ISFFT_Allocation.o \
  $(BUILD_DIR)/ISFFT_Planning.o $(BUILD_DIR)/ISFFT_Precompute.o
$(BUILD_DIR)/ISFFT_Redistribution.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Problem_State.o \
  $(BUILD_DIR)/ISFFT_MPI_Part.o $(BUILD_DIR)/ISFFT_Communication.o
$(BUILD_DIR)/ISFFT_Integrated_Solver.o: $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Plans.o \
  $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Workspace.o $(BUILD_DIR)/ISFFT_Halos_Exchange.o \
  $(BUILD_DIR)/ISFFT_Modal_Operations.o $(BUILD_DIR)/ISFFT_Communication.o $(BUILD_DIR)/ISFFT_Transform_Execution.o
$(BUILD_DIR)/ISFFT.o: $(BUILD_DIR)/ISFFT_Problem_State.o $(BUILD_DIR)/ISFFT_Initialization.o \
  $(BUILD_DIR)/ISFFT_Precision.o $(BUILD_DIR)/ISFFT_Integrated_Solver.o $(BUILD_DIR)/ISFFT_Transform_Execution.o \
  $(BUILD_DIR)/ISFFT_Redistribution.o

$(LIB_DIR)/$(LIB_NAME): $(ISFFT_OBJS) Makefile | $(LIB_DIR)
	@rm -f $@
	@$(AR) rcs $@ $(ISFFT_OBJS)

print-prefix:
	@printf 'PREFIX = %s\n' "$(PREFIX)"

install: library
	@mkdir -p "$(PREFIX)/lib" "$(PREFIX)/include"
	@cp $(LIB_DIR)/$(LIB_NAME) "$(PREFIX)/lib/"
	@cp $(INCLUDE_DIR)/*.mod "$(PREFIX)/include/"
	@printf 'ISFFT installed in %s\n' "$(PREFIX)"

clean:
	@rm -rf $(BUILD_DIR) $(INCLUDE_DIR) $(LIB_DIR)
