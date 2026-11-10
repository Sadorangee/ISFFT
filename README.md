# ISFFT

ISFFT is a self-contained, MPI-parallel FFT library with an integrated Poisson
solver for one-, two- and three-dimensional problems. It provides its own FFT,
DST and DCT kernels and uses MPI for parallel communication.

## Features

- Standalone real-to-complex and complex-to-complex FFT interfaces, with in-place
  and out-of-place variants.
- Poisson solvers with periodic, homogeneous Dirichlet, homogeneous Neumann and
  mixed boundary conditions.
- Spectral, explicit central-difference and compact central-difference schemes.
- Block, slab and pencil decompositions, with round-trip redistribution between
  the original decomposition and the internal transform decomposition.
- Arrays with or without halo cells.
- Estimate and Measure planning, with reusable execution workspaces.

## Project layout

```text
ISFFT_Github/
├── src/               ISFFT library sources
├── example/           CFD toy solver, input, MATLAB script and Makefile
├── Makefile           Library build and installation
├── CMakeLists.txt     CMake build and installation
├── README.md
└── LICENSE
```

Running `make` in the project root creates

```text
build/                 Library object files (.o)
include/               Fortran module files (.mod)
lib/libisfft.a         Static library
```

The library has 22 source files. The main groups are

- Basic state and storage in `Communicator`, `Precision`, `Problem_State`,
  `Discrete_Operators`, `Workspace` and `Plans`.
- MPI partitioning, halo exchange, transform kernels and modal operations.
- Communication, single-axis execution and complete transform execution.
- Allocation, precomputation and planning, coordinated by `Initialization`.
- Public redistribution and integrated Poisson interfaces, exported through
  the `ISFFT` module.

## Requirements and tested environments

Building requires a Fortran compiler, an MPI implementation with the `mpi`
Fortran module, `ar` and either GNU Make or CMake 3.16 or newer with a supported
build tool. The commands below use `mpif90` and `mpirun` from the same MPI
installation. MATLAB is optional for post-processing.

The Make build, static linking and installation workflow was checked on
macOS 26.6.2 (Apple Silicon) with GNU Fortran 15.2.0, Open MPI 5.0.9 and GNU Make
3.81, using `-O3`.
The CMake workflow was checked in the same environment with CMake 4.4.2,
the Unix Makefiles generator and Release mode.

The synchronized CPC library sources were also tested on Red Hat Enterprise
Linux 7.7 (x86_64) with GNU Fortran 4.8.5 and Intel MPI 2021.14, using `-O3`.

## Build and run with Make

From the project root

```sh
make
cd example
make
mpirun -n 4 ./Poisson.out
```

The root Makefile builds the library. The example Makefile compiles `menu.f90`,
which includes the other toy-solver sources, and links it with `libisfft.a`.
Build the library first and rerun the example's `make` after rebuilding it.

The supplied input uses a 128 by 128 by 128 periodic grid and a 2 by 2 by 1
process grid, so it runs with four MPI processes. The process count must equal
`npx0*npy0*npz0` in `example/Poisson_solver_UI.in`.

The executable reads the input from its current working directory and writes
`Poisson_1D_IO.dat`, `Poisson_2D_IO.dat` or `Poisson_3D_IO.dat`.

Compiler settings can be supplied on the command line

```sh
make clean
make FC=mpif90 FFLAGS="-O3"
make -C example clean
make -C example FC=mpif90 FFLAGS="-O3"
```

Clean and rebuild both parts when changing the compiler, MPI installation or
compiler flags. Fortran `.mod` files must match the compiler used by the caller.
The default active compiler flag is `-O3`.

## Build and run with CMake

From the project root, build the static library with

```sh
cmake -S . -B cmake-build -DCMAKE_Fortran_COMPILER=mpif90
cmake --build cmake-build --parallel 4
```

The default build type is `Release`. The static library is written to
`cmake-build/lib/libisfft.a` and the library modules to `cmake-build/include/`.
Object files and generated build rules are also kept inside `cmake-build/`.

Enable the toy solver with

```sh
cmake -S . -B cmake-build -DISFFT_BUILD_EXAMPLE=ON
cmake --build cmake-build --parallel 4
cd cmake-build/example
mpirun -n 4 ./Poisson.out
```

Configuration copies the input file and `Post_Poisson.m` from `example/` to
`cmake-build/example/`. Edit the input and post-processing settings in the
source `example/` directory, then rebuild to refresh these copies. Run the
executable and MATLAB script from `cmake-build/example/`.

Install the library and all library module files with

```sh
# From the project root
cmake --install cmake-build
# Or select another installation directory
cmake --install cmake-build --prefix "$HOME/.local/isfft"
```

The default installation directory is the project's `install/`. The installed
`lib/` and `include/` directories can be used with the external-linking commands
below. The example Makefile also accepts `ISFFT_PREFIX=../cmake-build` to link
against the CMake build directly.

Choose the compiler when creating the build directory. Use a new build directory
when changing the compiler or MPI installation.

## Installation and external linking

### Use the project build directly

After the root `make`, another solver can link against the project directories

```sh
mpif90 -O3 -I/path/to/ISFFT_Github/include \
    your_solver.f90 \
    -L/path/to/ISFFT_Github/lib -lisfft \
    -o solver.out
```

Place the library after the source or object files that use it. The supplied
example uses this mode by default.

### Install inside the project

```sh
make install
```

The default prefix is the project's `install/` directory. Installation copies
`libisfft.a` to `install/lib/` and all library `.mod` files to `install/include/`.
Use `make print-prefix` to display the destination.

The same example can use this installation

```sh
cd example
make clean
make ISFFT_PREFIX=../install
mpirun -n 4 ./Poisson.out
```

### Install elsewhere on the computer

Choose a directory with `PREFIX`

```sh
make install PREFIX="$HOME/.local/isfft"
```

A shared installation can use a prefix such as `/usr/local`, provided the user
has permission to write there. The directory structure is the same for every
prefix.

An external solver can use the installed library with

```sh
ISFFT_PREFIX="$HOME/.local/isfft"
mpif90 -O3 -I"$ISFFT_PREFIX/include" \
    your_solver.f90 \
    -L"$ISFFT_PREFIX/lib" -lisfft \
    -o solver.out
```

To build the supplied example against that prefix

```sh
cd example
make clean
make ISFFT_PREFIX="$HOME/.local/isfft"
```

The example can also be copied to another project and built with the same
`ISFFT_PREFIX` setting. Clean the example before switching library prefixes.

## Library interfaces

Applications access the public library interfaces with

```fortran
use ISFFT
```

| Interface | Purpose |
| --- | --- |
| `ISFFT_Initialize` | Prepare geometry, storage, precomputed quantities and execution strategies. |
| `ISFFT_Finalize` | Release the initialized problem and cached resources. |
| `ISFFT_Get_Local_Layout` | Query local physical sizes and global starting indices. |
| `ISFFT_Poisson_solver` | Solve a Poisson problem using the selected discretization. |
| `ISFFT_Redistribution_Forward`, `ISFFT_Redistribution_Inverse` | Transfer between original and internal decompositions. |
| `ISFFT_1D_FFT_*`, `ISFFT_2D_FFT_*`, `ISFFT_3D_FFT_*` | Execute public FFTs, including `r2c`, `c2c` and their `_inplace` variants. |
| `ISFFT_Get_FFT_Layout`, `ISFFT_Get_FFT_Local_Sizes` | Query the sizes and organization required by FFT arrays. |

Initialize MPI before using ISFFT. Call `ISFFT_Initialize` before execution and
`ISFFT_Finalize` before `MPI_Finalize`. Initialization is collective over the
selected communicator, which defaults to `MPI_COMM_WORLD`.

The optional initialization arguments `planning_level` and `prepare_scope`
control separate choices

- `ISFFT_PLAN_ESTIMATE` chooses strategies using preset rules.
- `ISFFT_PLAN_MEASURE` chooses strategies using initialization-time trials.
- `prepare_scope` selects which interfaces to prepare. The toy solver uses
  `ISFFT_USE_REDISTRIBUTION_POISSON` for redistribution and Poisson solving.
  `ISFFT_USE_ALL` prepares all supported interface families.

The defaults are `ISFFT_PLAN_MEASURE` and `ISFFT_USE_ALL`. Select the scope during
initialization and reuse the prepared strategies across repeated calls. All
ranks must agree on the planning settings.

The current initialization call and Poisson execution are shown in
[example/main_solver.f90](example/main_solver.f90). Use `ISFFT_REAL_KIND` for
real arrays passed to the library.

## Example input

Edit `example/Poisson_solver_UI.in` and rerun the executable. Recompilation is
needed only when source code changes.

| Input | Meaning |
| --- | --- |
| `nx_global`, `ny_global`, `nz_global` | Global grid sizes. |
| `Dimension` | Number of active directions, from 1 to 3. |
| `Boundary(1:3)` | Per-axis boundary types, with 0 periodic, 1 homogeneous Dirichlet and 2 homogeneous Neumann. |
| `npx0`, `npy0`, `npz0` | Process-grid dimensions in the original decomposition. |
| `LAP` | Halo width, with 0 for arrays without halo cells. |
| `SLx`, `SLy`, `SLz` | Positive domain lengths. |
| `Repeat_count` | Positive integer number of Poisson calls, with a default of 1. |
| `Decomposition` | Internal decomposition, with 0 block, 1 slab and 2 pencil. |
| `Scheme_type` | Discretization, as listed below. |

Set grid sizes and process counts to 1 along inactive directions. For 1D, use
one MPI process and `Decomposition=0`. For 2D, use `Decomposition=0` or `1`.
Pencil decomposition is available in 3D. The 2D slab path requires both active
grid sizes to be divisible by the MPI process count.

Periodic grids omit the repeated endpoint. Dirichlet and Neumann grids include
both endpoints. In 1D and 2D, the first periodic transform axis requires an even
grid size.

| `Scheme_type` | Discretization |
| --- | --- |
| 0 | Spectral |
| 1 | Second-order explicit central difference (E2) |
| 2 | Fourth-order explicit central difference (E4) |
| 3 | Fourth-order compact central difference (C4) |
| 4 | Sixth-order explicit central difference (E6) |
| 5 | Sixth-order compact central difference (C6) |
| 6 | Spectral-like compact central difference |

### Manufactured solutions

Each active boundary combination defines one analytical test, giving 3 cases
in 1D, 9 in 2D and 27 in 3D. The reference solution is

\[
\phi(\boldsymbol{x})=\prod_{d=1}^{D}g_d(x_d),
\qquad
f=-\left(\sum_{d=1}^{D}k_d^2\right)\phi.
\]

For a direction of length \(L\), the directional mode and wave number are

| Boundary | Mode \(g(x)\) | Wave number \(k\) |
| --- | --- | --- |
| Periodic | \(\cos(2\pi x/L)\) | \(2\pi/L\) |
| Homogeneous Dirichlet | \(\sin(\pi x/L)\) | \(\pi/L\) |
| Homogeneous Neumann | \(\cos(\pi x/L)\) | \(\pi/L\) |

The example solves \(\nabla^2\phi=f\). Its source uses the continuous Laplacian,
so finite-difference results include spatial discretization error when compared
with this analytical solution.

## Post-processing

In MATLAB, select `example/` as the current folder. Match `boundary`,
`dimension`, `grid_size` and `domain_size` at the top of `Post_Poisson.m` to the
input, then run

```matlab
Post_Poisson
```

The script plots a 1D profile with the analytical solution as a black line and
numerical values as markers. It prints the absolute `Linf` error, the
quadrature-weighted `L2` error and the maximum discrete residual `Residuals`.
For periodic/Neumann-only problems, it removes the arbitrary constant offset
before comparing solutions. Keep the input file matched to the output so that
residual evaluation uses the correct `Scheme_type`.

## Cleaning

From the root, `make clean` removes `build/`, `include/` and `lib/`. An installed
prefix remains available until it is removed separately.

From `example/`, `make clean` removes its objects, module files, executable and
numerical `.dat` output files. Library and example cleanup are separate.

For CMake, `cmake --build cmake-build --target clean` removes compiled files.
Remove `cmake-build/` to also remove its configuration, copied runtime files and
numerical output.

## Contact

Questions and bug reports can be sent to
[lizecheng@buaa.edu.cn](mailto:lizecheng@buaa.edu.cn).

## License

ISFFT is distributed under the MIT License. See [LICENSE](LICENSE).
