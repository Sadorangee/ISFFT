# ISFFT

ISFFT is a self-contained Fortran MPI library for FFT-based fast Poisson solvers in computational fluid dynamics (CFD) and high-performance computing (HPC) applications.

## Features

- MPI-parallel Poisson solvers for one-, two-, and three-dimensional problems.
- Periodic, homogeneous Dirichlet, homogeneous Neumann, and mixed boundary conditions.
- Pseudo-spectral and finite-difference formulations.
- Block, slab, and pencil decompositions, with reversible repartitioning between CFD and transform layouts.
- Interfaces that accept arrays with or without ghost cells.
- Public real-to-complex and complex-to-complex FFT interfaces.
- No dependency on an external FFT library.

## Directory Structure

```text
ISFFT/
├── examples/
│   └── Poisson/
│       ├── Poisson_solver_UI.in  # Example runtime configuration
│       ├── Poisson_source.f90    # Manufactured source terms
│       └── Post.m                # Unified smooth-case post-processing
├── src/                      # ISFFT library source files
├── Makefile                  # Library build and installation rules
└── LICENSE                   # MIT License
```

## Prerequisites

Before building ISFFT, make sure the following software is available:

- A modern Fortran compiler, such as `gfortran` or `ifort`.
- An MPI implementation, such as OpenMPI or MPICH.
- An MPI Fortran compiler wrapper, such as `mpif90`.
- GNU Make.
- MATLAB or GNU Octave is optional and is used only for example post-processing.

The commands below assume that `mpif90`, `mpirun`, and `make` are available in `PATH`.

## Building the Library

From the ISFFT root directory, run:

```bash
make
```

This command:

- compiles sources from `src/` into `build/`;
- writes Fortran module files to `include/`; and
- creates the static library `lib/libisfft.a`.

The files listed in `ISFFT_SRCS` in the root `Makefile` are in Fortran module-dependency order. Keep this order when compiling the sources manually.

To select another compiler or add compiler flags, override the Make variables on the command line, for example:

```bash
make FC=mpif90 FFLAGS="-O3 -g"
```

## Installation

### System-Style Install

Install to a standard system location (requires root/admin privileges):

```bash
make install
```

By default, `PREFIX=/usr/local`. This installs:
- `lib/libisfft.a` → `/usr/local/lib/`
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
- `lib/libisfft.a` → `$(PREFIX)/lib/`
- Module files → `$(PREFIX)/include/`

## Linking ISFFT to Another Solver

### Using Repository Build Artifacts

After running `make` in the ISFFT root directory, compile an external solver with:

```bash
mpif90 -I/path/to/ISFFT/include \
       your_solver.f90 \
       -L/path/to/ISFFT/lib -lisfft \
       -o solver.out
```

The example Makefile uses `../../include` and `../../lib` by default:

```bash
cd examples/Poisson
make
```

### Using an Installed Prefix

For a library installed under a custom prefix:

```bash
export ISFFT_PREFIX=/path/to/ISFFT/install
mpif90 -I$ISFFT_PREFIX/include \
       your_solver.f90 \
       -L$ISFFT_PREFIX/lib -lisfft \
       -o solver.out
```

The example Makefile switches to the installed-prefix configuration when `ISFFT_PREFIX` is supplied:

```bash
cd examples/Poisson
make ISFFT_PREFIX=/path/to/ISFFT/install
```

No Makefile editing is required.

## Building and Running the Example

Build the library first, then build the example:

```bash
# From the ISFFT root directory
make

cd examples/Poisson
make
mpirun -n 8 ./Poisson.out
```

The number passed to `mpirun -n` must equal:

```text
npx0 × npy0 × npz0
```

The supplied `Poisson_solver_UI.in` uses `2 × 2 × 2`, so its default run requires eight MPI processes.

The executable reads `Poisson_solver_UI.in` from the current working directory and writes one of:

- `Poisson_1D_IO.dat`
- `Poisson_2D_IO.dat`
- `Poisson_3D_IO.dat`

## Example Configuration

The example is controlled by `examples/Poisson/Poisson_solver_UI.in`.

### Grid, dimension, and boundary conditions

```text
nx_global, ny_global, nz_global, Dimension, Boundary(1:3)
```

| Parameter | Meaning |
| --- | --- |
| `nx_global`, `ny_global`, `nz_global` | Global grid sizes. |
| `Dimension` | Active dimension: `1`, `2`, or `3`. |
| `Boundary(d)=0` | Periodic boundary condition in direction `d`. |
| `Boundary(d)=1` | Homogeneous Dirichlet boundary condition. |
| `Boundary(d)=2` | Homogeneous Neumann boundary condition. |

Dimension-specific settings:

- `Dimension=1` requires `ny_global=nz_global=1`, `npy0=npz0=1`, one MPI process, and `If_decomp=0`.
- `Dimension=2` requires `nz_global=1`, `npz0=1`, and `If_decomp=0` or `1`.
- `Dimension=3` supports `If_decomp=0`, `1`, or `2`.
- For inactive directions, the corresponding boundary value is ignored.
- The first periodic transform direction currently requires an even number of grid points.

Periodic grids do not duplicate the endpoint. Dirichlet and Neumann grids include both boundary points.

### Process grid and ghost cells

```text
npx0, npy0, npz0, LAP
```

| Parameter | Meaning |
| --- | --- |
| `npx0`, `npy0`, `npz0` | MPI process-grid dimensions in the physical layout. Their product must equal the MPI process count. |
| `LAP` | Number of ghost-cell layers around each local array. Use `0` for arrays without ghost cells. |

Each process-grid dimension must not exceed the corresponding active global grid size.

### Computational domain

```text
SLx, SLy, SLz
```

These values specify the physical domain lengths. Inactive directions still require positive placeholder values.

### Time-control values in the example

```text
CFL, T
```

The example reads these values and computes a nominal time step. For the current solver demonstration, `main_solver.f90` explicitly sets the number of Poisson calls to one. These parameters do not affect the library interface itself.

### Decomposition type

| `If_decomp` | Layout | Availability |
| ---: | --- | --- |
| `0` | Block | 1D, 2D, and 3D |
| `1` | Slab | 2D and 3D |
| `2` | Pencil | 3D only |

When a caller's physical process grid differs from the transform layout, ISFFT performs reversible repartitioning internally for the Poisson solver.

### Discretization scheme

| `If_scheme` | Scheme |
| ---: | --- |
| `0` | Pseudo-spectral |
| `1` | Second-order central difference |
| `2` | Fourth-order central difference |
| `3` | Fourth-order compact central difference |
| `4` | Sixth-order central difference |
| `5` | Sixth-order compact central difference |
| `6` | Spectral-like compact central difference |

### Manufactured-test selectors

```text
If_Periodic_test, If_Dirichlet_test, If_Neumann_test
```

These selectors are used only for three-dimensional cases with the same boundary type in all three directions:

- `[0,0,0]` uses `If_Periodic_test`.
- `[1,1,1]` uses `If_Dirichlet_test`.
- `[2,2,2]` uses `If_Neumann_test`.
- 1D, 2D, and 3D mixed-boundary cases use fixed manufactured solutions, so the three selectors are ignored.

## Manufactured Test Cases

The example solves

```math
\nabla^2 \phi = f.
```

The source term `f` is generated analytically from the reference solution shown below.

### 1D and 2D cases

For 1D, the fixed reference solution is:

| Boundary | Reference solution |
| --- | --- |
| P | `sin(2πx/Lx)` |
| D | `sin(πx/Lx)` |
| N | `cos(πx/Lx)` |

For 2D, all nine combinations of P, D, and N are supported:

| Boundary | Reference solution |
| --- | --- |
| PP | `sin(2πx/Lx) + cos(4πy/Ly)` |
| PD | `cos(2πx/Lx) sin(πy/Ly)` |
| PN | `cos(2πx/Lx) cos(πy/Ly)` |
| DP | `sin(πx/Lx) cos(2πy/Ly)` |
| DD | `sin(πx/Lx) sin(2πy/Ly)` |
| DN | `sin(πx/Lx) cos(πy/Ly)` |
| NP | `cos(πx/Lx) cos(2πy/Ly)` |
| ND | `cos(πx/Lx) sin(πy/Ly)` |
| NN | `cos(πx/Lx) cos(2πy/Ly)` |

### 3D periodic cases

| `If_Periodic_test` | Reference solution |
| ---: | --- |
| `0` | `sin(2πz/Lz) [sin(2πx/Lx) + cos(4πy/Ly)]` |
| `1` | `x²(Lx-x)² y²(Ly-y)² z²(Lz-z)²`, with its mean removed |

### 3D Dirichlet cases

| `If_Dirichlet_test` | Reference or source |
| ---: | --- |
| `0` | `sin(πx/Lx) sin(2πy/Ly) sin(3πz/Lz)` |
| `1` | `sin(πx) sin(πy) sin(πz)` |
| `2` | `x(Lx-x) y(Ly-y) z(Lz-z)` |
| `3` | `x²(Lx-x)² y²(Ly-y)² z²(Lz-z)²` |
| `4` | `x³(Lx-x)³ y³(Ly-y)³ z³(Lz-z)³` |
| `5` | Gaussian multiplied by `sin²` boundary factors in each direction |
| `6` | Product of smooth tanh windows; it is near zero, but not exactly zero, on the boundary |
| `7` | Discrete point source under homogeneous Dirichlet boundary conditions |

### 3D Neumann cases

| `If_Neumann_test` | Reference or source |
| ---: | --- |
| `0` | `cos(πx) cos(2πy) cos(3πz)` |
| `1` | `x²(Lx-x)² y²(Ly-y)² z²(Lz-z)²` |
| `2` | Discrete point source with a zero-mean compatibility correction |

### 3D mixed-boundary cases

For any 3D combination that is not PPP, DDD, or NNN, the reference solution is the product of one factor per direction:

```text
Periodic:   cos(2πx/L)
Dirichlet:  sin(πx/L)
Neumann:    cos(πx/L)
```

The coordinate and domain length corresponding to each direction are used in its factor.

## Post-Processing

For all smooth manufactured cases, open `examples/Poisson/Post.m` and set the values near the top to match `Poisson_solver_UI.in`:

```matlab
Dimension = 3;
Boundary = [0, 0, 0];
If_case = 0;
Grid_Size = [128, 128, 128];
Domain_Size = [4.0, 4.0, 4.0];
```

Then run `Post.m` from the example directory. It:

- loads the dimension-appropriate output file;
- reconstructs the exact manufactured solution;
- removes the arbitrary constant offset for problems without a Dirichlet boundary;
- reports the relative L2 error, `||phi_ISFFT-phi_exact||₂ / ||phi_exact||₂`;
- reports the absolute L-infinity error; and
- plots a representative solution profile against the exact solution.

The singular point-source cases require different diagnostics and retain their specialized scripts:

- Dirichlet case 7: `Post_Dirichlet_cases.m`
- Neumann case 2: `Post_Neumann_cases.m`

Set the corresponding `If_case` and grid/domain values near the top of those scripts before running them.

## Cleaning

Clean the library build artifacts from the root directory:

```bash
make clean
```

Clean example object files, output data, and the executable separately:

```bash
cd examples/Poisson
make clean
```

## Contact

For questions and bug reports, contact [lizecheng@buaa.edu.cn](mailto:lizecheng@buaa.edu.cn).

## License

ISFFT is distributed under the MIT License. See [LICENSE](LICENSE) for details.
