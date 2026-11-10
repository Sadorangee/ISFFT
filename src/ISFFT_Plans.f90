!------------------------------------------------------------------------------
! ISFFT_Plans
! Selected execution strategies and their stage descriptions.
!------------------------------------------------------------------------------
module ISFFT_Plans
   use mpi, only: MPI_COMM_NULL
   use ISFFT_Precision, only: ISFFT_INT64_KIND
   use ISFFT_Problem_State, only: Layout
   implicit none
   private
   public :: Axis_Plan, axis_plan_cache, block_cache, Block_Choice, default_chunk_elements
   public :: default_stream_elements, Kernel_Plan, poisson_plan, prepared_kernel_plans, public_export
   public :: public_import, public_plans, redistribution_forward_plan, redistribution_inverse_plan, Transfer_Plan
   public :: Transform_Plan

   integer(ISFFT_INT64_KIND), parameter :: default_chunk_elements = 67108864_ISFFT_INT64_KIND
   integer(ISFFT_INT64_KIND), parameter :: default_stream_elements = 1048576_ISFFT_INT64_KIND

   type :: Transfer_Plan
      integer :: resource_id = 0
      type(Layout) :: source, destination
      !! Counts and packed displacements are in field elements.
      integer(ISFFT_INT64_KIND), allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      !! Box offsets are zero-based and use physical-axis order.
      integer, allocatable :: send_box_offsets(:, :), recv_box_offsets(:, :), send_box_sizes(:, :), recv_box_sizes(:, :)
      integer, allocatable :: send_types(:), recv_types(:), send_type_counts(:), recv_type_counts(:), &
         send_byte_displs(:), recv_byte_displs(:)
      integer :: gather_comm = MPI_COMM_NULL, exchange_comm = MPI_COMM_NULL
      integer(ISFFT_INT64_KIND) :: gather_count = 0, max_peer = 0, chunk_limit = default_chunk_elements
      !! chunked is a mandatory count-safety constraint, not a tuned choice.
      logical :: chunked = .false., streamed = .false.
      integer(ISFFT_INT64_KIND) :: stream_width = default_stream_elements
      integer, allocatable :: peers(:)
      integer(ISFFT_INT64_KIND), allocatable :: gather_counts(:), gather_displs(:)
      logical :: packed_p2p_active = .false.
      logical :: direct = .false., sparse = .false., local_only = .false., ready = .false.
      !! Local exchange-communicator ranks; global datatype indices are peers(list).
      !! Includes self, using disjoint input/output buffers. Requests complete per call.
      integer, allocatable :: sparse_send(:), sparse_recv(:), requests(:)
   end type

   !! Numerical, field and transform plan records
   !! transform_kind selects FFT-R2C (1), FFT-C2C (2), DST-I (3) or DCT-I (4).
   !! kernel_method selects scalar/batched kernel implementations (kernel-specific codes).
   !! inner_batch_width is the butterfly/real-transform microbatch width.
   !! outer_batch_width is the independently tuned number of lines submitted per call.
   !! Axis_Plan%communication_batch_width independently controls MPI block tiles.
   type :: Kernel_Plan
      integer :: transform_length = 0, transform_kind = 0, line_count = 0
      integer :: kernel_method = 1, inner_batch_width = 1, outer_batch_width = 1
      integer :: local_shape(3) = 0, memory_axis = 1
      !! Access selection is independent of the positive widths used by either path.
      logical :: direct_access = .false.
      integer :: direct_batch_width = 1 !! Lines per direct strided FFT microbatch.
      logical :: complex_input = .false.
   end type

   type :: Axis_Plan
      integer :: physical_axis = 0, transform_kind = 0, transform_length = 0, communication_batch_width = 1
      logical :: tiled_block = .false.
      type(Kernel_Plan) :: tail_kernel
      type(Kernel_Plan) :: kernel
      type(Layout) :: input, output
      type(Transfer_Plan) :: before, after, reverse_before, reverse_after, pack, unpack, reverse_pack, reverse_unpack
   end type

   type :: Transform_Plan
      integer :: dimensions = 0, half_axis = 0, generation = -1
      type(Layout) :: physical, spectral
      type(Axis_Plan) :: axes(3)
      logical :: ready = .false.
      integer :: layout_choice = 0
   end type

   type :: Block_Choice
      type(Axis_Plan) :: plan
      integer :: method, width, line_count
   end type

   type(Kernel_Plan), allocatable, save :: prepared_kernel_plans(:), axis_plan_cache(:)
   type(Transform_Plan), save :: poisson_plan, public_plans(2)
   type(Transfer_Plan), save :: redistribution_forward_plan, redistribution_inverse_plan
   type(Transfer_Plan), save :: public_export(2, 2), public_import(2, 2)
   type(Block_Choice), allocatable, save :: block_cache(:)

contains


end module ISFFT_Plans
