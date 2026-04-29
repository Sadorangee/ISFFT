# Example Solver Makefile for ISFFT
# This Makefile compiles the example solver (menu.f90) and links against the ISFFT library.
# It assumes ISFFT has already been built or installed.

# Compiler settings
FC = mpif90
#FC = ifort
FFLAGS = -O3 -g

# Executable name
EXE = 3DPeriodicPoisson.out

# ----------------------------------------------------------------------
# ISFFT library configuration
#
# Default mode:
#   use ISFFT directly from repository build artifacts
#
# Installed-prefix mode:
#   make ISFFT_PREFIX=/path/to/ISFFT/install
# ----------------------------------------------------------------------

ifeq ($(origin ISFFT_PREFIX), undefined)
ISFFT_INCLUDE = ../../include
ISFFT_LIB = ../../lib
ISFFT_MODE = repository build artifacts
else
ISFFT_INCLUDE = $(ISFFT_PREFIX)/include
ISFFT_LIB = $(ISFFT_PREFIX)/lib
ISFFT_MODE = installed prefix: $(ISFFT_PREFIX)
endif

ISFFT_STATIC_LIB = $(ISFFT_LIB)/libisfft.a

CPPFLAGS = -I$(ISFFT_INCLUDE)
LDFLAGS  = -L$(ISFFT_LIB)
LIBS     = -lisfft

.PHONY: all clean print-config check-isfft

all: check-isfft $(EXE)

print-config:
	@echo "ISFFT mode    = $(ISFFT_MODE)"
	@echo "ISFFT_INCLUDE = $(ISFFT_INCLUDE)"
	@echo "ISFFT_LIB     = $(ISFFT_LIB)"
	@echo "EXE            = $(EXE)"


# Friendly checks before compilation
check-isfft:
	@if [ ! -d "$(ISFFT_INCLUDE)" ]; then \
		echo "Error: ISFFT include directory not found: $(ISFFT_INCLUDE)"; \
		echo "Hint: build ISFFT first with 'make' in the ISFFT root directory,"; \
		echo "      or use 'make ISFFT_PREFIX=/path/to/install' if you want installed mode."; \
		exit 1; \
	fi
	@if [ ! -f "$(ISFFT_STATIC_LIB)" ]; then \
		echo "Error: ISFFT static library not found: $(ISFFT_STATIC_LIB)"; \
		echo "Hint: run 'make' in the ISFFT root directory to generate lib/libisfft.a,"; \
		echo "      or install ISFFT and pass ISFFT_PREFIX=/path/to/install."; \
		exit 1; \
	fi
	@if ! ls "$(ISFFT_INCLUDE)"/*.mod >/dev/null 2>&1; then \
		echo "Error: no .mod files found in $(ISFFT_INCLUDE)"; \
		echo "Hint: ISFFT must export Fortran module files before this example can compile."; \
		echo "      Run 'make' in the ISFFT root directory first."; \
		exit 1; \
	fi

$(EXE): menu.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod *.dat *.out $(EXE)
