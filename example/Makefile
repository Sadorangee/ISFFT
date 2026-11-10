FC = mpif90
FFLAGS = -O3

# Use the repository build by default, or select an installed prefix.
ISFFT_PREFIX ?= ..
ISFFT_INCLUDE = $(ISFFT_PREFIX)/include
ISFFT_STATIC_LIB = $(ISFFT_PREFIX)/lib/libisfft.a
ISFFT_MODS = $(wildcard $(ISFFT_INCLUDE)/*.mod)
SOLVER_SRCS = menu.f90 parameters.f90 part_mpi.f90 main_solver.f90 Poisson_source.f90 IO_Poisson.f90
EXE = Poisson.out

.NOTPARALLEL:
.DELETE_ON_ERROR:

.PHONY: all clean check-library
all: check-library $(EXE)

check-library:
	@test -f "$(ISFFT_STATIC_LIB)" -a -f "$(ISFFT_INCLUDE)/isfft.mod" || { \
	  printf 'Build ISFFT in the project root, or set ISFFT_PREFIX to its installation.\n' >&2; exit 1; }
ifeq ($(abspath $(ISFFT_PREFIX)),$(abspath ..))
	@$(MAKE) --no-print-directory -C .. -q lib/libisfft.a || { \
	  printf 'Rebuild ISFFT by running make in the project root.\n' >&2; exit 1; }
endif

# menu.f90 includes all external solver sources.
menu.o: $(SOLVER_SRCS) $(ISFFT_MODS) Makefile | check-library
	@$(FC) $(FFLAGS) -I"$(ISFFT_INCLUDE)" -J. -c menu.f90 -o $@

$(EXE): menu.o $(ISFFT_STATIC_LIB) Makefile | check-library
	@$(FC) $(FFLAGS) menu.o "$(ISFFT_STATIC_LIB)" -o $@
	@printf 'CFD solver built successfully.\n'

clean:
	@rm -f *.o *.mod *.out *.dat $(EXE)
