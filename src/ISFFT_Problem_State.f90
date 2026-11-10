!------------------------------------------------------------------------------
! ISFFT_Problem_State
! Problem geometry, data descriptions and initialization settings.
!------------------------------------------------------------------------------
module ISFFT_Problem_State
   use mpi
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   implicit none
   private
   public :: boundary_types, decomposition_type, dimensions, discretization_scheme, domain_length_x
   public :: domain_length_y, domain_length_z, fft_output_modes, halfspectrum_x, halfspectrum_y
   public :: halfspectrum_z, halo_width, i_nn, i_offset, ID_XM1
   public :: ID_XP1, ID_YM1, ID_YP1, ID_ZM1, ID_ZP1
   public :: internal_layout, ISFFT_Barrier_level, ISFFT_FFT_BOTH, ISFFT_FFT_INTERNAL, ISFFT_FFT_ORIGINAL
   public :: ISFFT_PLAN_ESTIMATE, ISFFT_plan_generation, ISFFT_PLAN_MEASURE, ISFFT_USE_ALL, ISFFT_USE_C2CFFT_INPLACE
   public :: ISFFT_USE_C2CFFT_OUTOFPLACE, ISFFT_USE_POISSON, ISFFT_USE_R2CFFT_INPLACE, &
      ISFFT_USE_R2CFFT_OUTOFPLACE, ISFFT_USE_REDISTRIBUTION
   public :: ISFFT_USE_REDISTRIBUTION_C2CFFT_INPLACE, ISFFT_USE_REDISTRIBUTION_C2CFFT_OUTOFPLACE, &
      ISFFT_USE_REDISTRIBUTION_POISSON, ISFFT_USE_REDISTRIBUTION_R2CFFT_INPLACE, &
      ISFFT_USE_REDISTRIBUTION_R2CFFT_OUTOFPLACE
   public :: j_nn, j_offset, k_nn, k_offset, Layout
   public :: layout_policy, MPI_COMM_X, MPI_COMM_XY, MPI_COMM_XZ, MPI_COMM_Y
   public :: MPI_COMM_YZ, MPI_COMM_Z, my_id, np_size, npx
   public :: npx0, npy, npy0, npz, npz0
   public :: nx, nx_global, ny, ny_global, nz
   public :: nz_global, original_layout, planning_level, prepare_scope, public_fft_modes
   public :: public_layout, topology_ready, transport_policy

   !! Geometry, boundary conditions and problem generation
   integer, save :: decomposition_type, discretization_scheme, dimensions
   integer, save :: ISFFT_plan_generation = 0
   integer, save :: nx_global, ny_global, nz_global, nx, ny, nz, halo_width
   integer, save :: boundary_types(3)
   real(kind=ISFFT_REAL_KIND), save :: domain_length_x, domain_length_y, domain_length_z

   !! Active process coordinates, neighbors and communicator handles
   !! Switch ownership through the layout activation routines.
   integer, save :: my_id, npx, npy, npz, npx0, npy0, npz0, ID_XP1, ID_XM1, ID_YP1, ID_YM1, ID_ZP1, ID_ZM1, &
      MPI_COMM_X, MPI_COMM_Y, MPI_COMM_Z, MPI_COMM_XY, MPI_COMM_XZ, MPI_COMM_YZ, ISFFT_Barrier_level, np_size

   !! Process indices start at zero. Offsets store one-based global starts.
   !! Keep table capacity consistent with the guard in ISFFT_part.
   integer, save :: i_offset(0:2048), j_offset(0:2048), k_offset(0:2048), i_nn(0:2048), j_nn(0:2048), k_nn(0:2048)

   !! Initialization settings and prepared interface modes
   !! Public planning flags and usage bit assignments
   !! Keep bit values aligned with initialization validation and public API exports.
   !! ESTIMATE: no performance trials; MEASURE: sampled candidate screening.
   !! Usage bits: redistribution, Poisson, R2C internal/original, C2C internal/original.
   integer, parameter :: ISFFT_PLAN_ESTIMATE = 0, ISFFT_PLAN_MEASURE = 1
   integer, parameter :: ISFFT_USE_REDISTRIBUTION = 1
   integer, parameter :: ISFFT_USE_POISSON = 2
   integer, parameter :: ISFFT_USE_REDISTRIBUTION_POISSON = 3
   integer, parameter :: ISFFT_USE_R2CFFT_OUTOFPLACE = 4
   integer, parameter :: ISFFT_USE_R2CFFT_INPLACE = 8
   integer, parameter :: ISFFT_USE_C2CFFT_OUTOFPLACE = 16
   integer, parameter :: ISFFT_USE_C2CFFT_INPLACE = 32
   integer, parameter :: ISFFT_USE_REDISTRIBUTION_R2CFFT_OUTOFPLACE = 5
   integer, parameter :: ISFFT_USE_REDISTRIBUTION_R2CFFT_INPLACE = 9
   integer, parameter :: ISFFT_USE_REDISTRIBUTION_C2CFFT_OUTOFPLACE = 17
   integer, parameter :: ISFFT_USE_REDISTRIBUTION_C2CFFT_INPLACE = 33
   integer, parameter :: ISFFT_USE_ALL = 63
   integer, save :: public_fft_modes(2) = 0
   integer, parameter :: ISFFT_FFT_INTERNAL = 1, ISFFT_FFT_ORIGINAL = 2, ISFFT_FFT_BOTH = 3
   integer, save :: fft_output_modes = ISFFT_FFT_BOTH
   integer, save :: planning_level = ISFFT_PLAN_MEASURE, prepare_scope = ISFFT_USE_ALL
   logical, save :: topology_ready = .false.

   !! Data decomposition and memory organization
   !! Geometry and boxes use physical-axis order.
   !! axis_order maps memory axes to physical axes.
   type :: Layout
      integer :: global_shape(3) = 1, process_grid(3) = 1, axis_order(3) = [1, 2, 3]
      integer :: replicated_axis = 0, partition_shape(3) = 0
      logical :: complex_data = .false.
      !! rank_boxes stores sizes in rows 1-3 and zero-based offsets in rows 4-6.
      integer, allocatable :: rank_boxes(:, :)
   end type

   type(Layout), save :: original_layout, internal_layout, public_layout(2, 2)
   integer, save :: halfspectrum_x = 0, halfspectrum_y = 0, halfspectrum_z = 0
   integer, save :: layout_policy = -1, transport_policy = -1

contains


end module ISFFT_Problem_State
