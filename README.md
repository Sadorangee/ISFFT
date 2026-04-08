# SPLASH

SPLASH is a Fortran MPI library for Poisson equation solvers in computational fluid dynamics (CFD) and high-performance computing (HPC) applications.

## Features

SPLASH is designed as a plug-and-play FFT-based fast Poisson solver for block-decomposed CFD workflows. Its main features include:

- **MPI-parallel 3D Poisson solver** for large-scale CFD and HPC applications.
- **FFT-based fast solution strategy** with support for both pseudo-spectral and finite-difference formulations under periodic boundary conditions.
- **Reversible repartitioning framework** that bridges the mismatch between native block decomposition in CFD codes and FFT-compatible slab/pencil layouts.
- **Minimal integration into existing solvers**, allowing users to call the Poisson solver without reorganizing the host code's overall workflow.
- **Ghost-cell-aware interface**, so padded CFD arrays can be passed directly to the solver.
- **Self-contained software design**, avoiding external FFT library dependencies and simplifying deployment.

## Directory Structure

```text
SPLASH/
├── examples/       # Example solvers using the SPLASH library
│   └── 3D_Periodic_Poisson/ # 3D Poisson solver example
├── src/            # SPLASH library source files
├── Makefile        # Main Makefile to build the library
└── LICENSE         # MIT License file
```

## Prerequisites

Before building SPLASH, make sure the following software is available in your environment:

- **Fortran compiler** with modern Fortran support, such as `gfortran` or `ifort`.
- **MPI implementation** providing the MPI compiler wrapper and runtime, such as OpenMPI or MPICH.
- **GNU Make** for building the library and examples.

> The example commands in this README assume an MPI-enabled Fortran compiler command such as `mpif90` is available in your `PATH`.

## Building the Library

The default `make` command builds only the SPLASH static library:

```bash
make
```

This compiles source files from `src/`, places object files in `build/`, module files in `include/`, and creates the static library `lib/libsplash.a`.

**Important for Source-Level Development**: If you are modifying the SPLASH root source code and compiling it manually (not using the provided Makefile), Fortran module dependencies require files to be compiled in a specific order. Please refer to the `SPLASH_SRCS` variable in the `Makefile` for the exact dependency order.

## Installation

### System-Style Install

Install to a standard system location (requires root/admin privileges):

```bash
make install
```

By default, `PREFIX=/usr/local`. This installs:
- `lib/libsplash.a` → `/usr/local/lib/`
- Module files → `/usr/local/include/`

### Local Project Install

Install to a custom directory (no special privileges needed):

```bash
make install PREFIX=/path/to/your/project/install
```

For example, to install into a local `install` folder in the current directory:

```bash
make install PREFIX=$(PWD)/install
```

This creates:
- `lib/libsplash.a` → `$(PREFIX)/lib/`
- Module files → `$(PREFIX)/include/`

## Using SPLASH in External Solvers

External solvers can link against SPLASH in two ways:

### Mode A: Link Against Repository Build Artifacts

If SPLASH is built but not installed, use relative paths:

```bash
mpif90 -I/path/to/SPLASH/include -L/path/to/SPLASH/lib -lsplash your_solver.f90 -o solver.out
```

For the example solver in `examples/3D_Periodic_Poisson/`:

```bash
cd examples/3D_Periodic_Poisson
make  # Uses ../../include and ../../lib by default
```

### Mode B: Link Against Installed SPLASH

If SPLASH is installed to a prefix, use the prefix paths:

```bash
export SPLASH_PREFIX=/usr/local  # or your custom PREFIX
mpif90 -I$SPLASH_PREFIX/include -L$SPLASH_PREFIX/lib -lsplash your_solver.f90 -o solver.out
```

For the example solver using installed SPLASH:

```bash
cd examples/3D_Periodic_Poisson
# Edit Makefile to uncomment Mode B and set SPLASH_PREFIX
make
```

## Building and Running the Example Solver

The example solver demonstrates Poisson 3D equation solving:

```bash
# Build SPLASH library first
make

# Build and run the example
cd examples/3D_Periodic_Poisson
make
mpirun -n 4 ./3DPeriodicPoisson.out  # Runs with 4 processes (adjust -n as needed)
```

### Customizing the Example

- **Source Term**: Defined in `Poisson_3D_source.f90`. Modify this file to change the source term (requires recompilation).
- **Parameters**: Modify `Poisson_solver_UI.in` to test different configurations (e.g., grid size) *without* recompiling.
- **Post-processing**: Use `Post.m` in MATLAB or Octave to visualize or process the results.

## Cleaning

Clean intermediate and built artifacts. You need to clean the main library and the example separately:

```bash
# In the root directory (cleans library build/, include/, and lib/):
make clean

# In the example directory (cleans example object files and executable):
cd examples/3D_Periodic_Poisson
make clean
```

## Contact

If you have any questions or bug reports, please contact [lizecheng@buaa.edu.cn](mailto:lizecheng@buaa.edu.cn).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.