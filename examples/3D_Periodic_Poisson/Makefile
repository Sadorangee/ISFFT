# Example Solver Makefile for SPLASH
# This Makefile compiles the example solver (menu.f90) and links against the SPLASH library.
# It assumes SPLASH has already been built or installed.

# Compiler settings
FC = mpif90
#FC = ifort
FFLAGS = -O3 -g

# Executable name
EXE = 3DPeriodicPoisson.out

# ----------------------------------------------------------------------
# SPLASH library configuration
#
# Default mode:
#   use SPLASH directly from repository build artifacts
#
# Installed-prefix mode:
#   make SPLASH_PREFIX=/path/to/SPLASH/install
# ----------------------------------------------------------------------

ifeq ($(origin SPLASH_PREFIX), undefined)
SPLASH_INCLUDE = ../../include
SPLASH_LIB = ../../lib
SPLASH_MODE = repository build artifacts
else
SPLASH_INCLUDE = $(SPLASH_PREFIX)/include
SPLASH_LIB = $(SPLASH_PREFIX)/lib
SPLASH_MODE = installed prefix: $(SPLASH_PREFIX)
endif

SPLASH_STATIC_LIB = $(SPLASH_LIB)/libsplash.a

CPPFLAGS = -I$(SPLASH_INCLUDE)
LDFLAGS  = -L$(SPLASH_LIB)
LIBS     = -lsplash

.PHONY: all clean print-config check-splash

all: check-splash $(EXE)

print-config:
	@echo "SPLASH mode    = $(SPLASH_MODE)"
	@echo "SPLASH_INCLUDE = $(SPLASH_INCLUDE)"
	@echo "SPLASH_LIB     = $(SPLASH_LIB)"
	@echo "EXE            = $(EXE)"


# Friendly checks before compilation
check-splash:
	@if [ ! -d "$(SPLASH_INCLUDE)" ]; then \
		echo "Error: SPLASH include directory not found: $(SPLASH_INCLUDE)"; \
		echo "Hint: build SPLASH first with 'make' in the SPLASH root directory,"; \
		echo "      or use 'make SPLASH_PREFIX=/path/to/install' if you want installed mode."; \
		exit 1; \
	fi
	@if [ ! -f "$(SPLASH_STATIC_LIB)" ]; then \
		echo "Error: SPLASH static library not found: $(SPLASH_STATIC_LIB)"; \
		echo "Hint: run 'make' in the SPLASH root directory to generate lib/libsplash.a,"; \
		echo "      or install SPLASH and pass SPLASH_PREFIX=/path/to/install."; \
		exit 1; \
	fi
	@if ! ls "$(SPLASH_INCLUDE)"/*.mod >/dev/null 2>&1; then \
		echo "Error: no .mod files found in $(SPLASH_INCLUDE)"; \
		echo "Hint: SPLASH must export Fortran module files before this example can compile."; \
		echo "      Run 'make' in the SPLASH root directory first."; \
		exit 1; \
	fi

$(EXE): menu.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod *.dat $(EXE)