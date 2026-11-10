!------------------------------------------------------------------------------
! ISFFT_Initialization
! Public problem initialization, finalization and state queries.
!------------------------------------------------------------------------------
module ISFFT_Initialization
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm, ISFFT_Release_Communicator, ISFFT_Set_Communicator
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   use ISFFT_Plans, only: axis_plan_cache, block_cache
   use ISFFT_Problem_State, only: domain_length_x, domain_length_y, domain_length_z, fft_output_modes, &
      halfspectrum_x, halfspectrum_y, halfspectrum_z, i_offset, &
      ISFFT_FFT_BOTH, ISFFT_plan_generation, ISFFT_PLAN_MEASURE, ISFFT_USE_ALL, &
      j_offset, k_offset, layout_policy, my_id, &
      np_size, npx, npx0, npy, &
      npy0, npz, npz0, nx, &
      nx_global, ny, ny_global, nz, &
      nz_global, public_fft_modes, topology_ready, transport_policy, &
      configured_planning_level => planning_level, configured_prepare_scope => prepare_scope
   use ISFFT_Workspace, only: Release_Execution_Scratch
   use ISFFT_MPI_Part, only: ISFFT_Get_Internal_Process_Grid, ISFFT_part, ISFFT_Part_Finalize
   use ISFFT_Communication, only: Release_Transport_Cache
   use ISFFT_Allocation, only: Finalize_Plans, Prepare_Execution_Workspaces
   use ISFFT_Planning, only: Ensure_Public_Plan, Initialize_Plans, Policy_Configure
   use ISFFT_Precompute, only: FFT_Finalize
   implicit none
   private
   public :: ISFFT_Finalize, ISFFT_Get_Local_Layout, ISFFT_Initialize, Set_Plan_Options

contains

   !! Collective initialization on ISFFT_Comm. Optional planning settings must agree on all ranks.
   !! Omitted settings preserve the callable API: MEASURE prepares Poisson, C2C and R2C plans.
   !! Select prepare_scope explicitly to avoid preparing unused public interfaces. Calls to an
   !! unprepared solver/transform fail; execution never triggers lazy planning.
   subroutine ISFFT_Initialize( &
      dimensions, global_grid_size, domain_lengths, &
      boundary_types, discretization_scheme, &
      process_grid_size, decomposition_type, halo_width, rank, &
      planning_level, prepare_scope, comm, fft_output)
      use ISFFT_Problem_State, only: &
         configured_decomposition_type => decomposition_type, &
         configured_discretization_scheme => discretization_scheme, &
         configured_dimensions => dimensions, &
         configured_boundary_types => boundary_types, &
         configured_halo_width => halo_width
      integer, intent(in) :: dimensions, global_grid_size(3)
      real(kind=ISFFT_REAL_KIND), intent(in) :: domain_lengths(3)
      integer, intent(in) :: boundary_types(3), discretization_scheme
      integer, intent(in) :: process_grid_size(3), decomposition_type, halo_width, rank
      integer, intent(in), optional :: planning_level, prepare_scope
      integer, intent(in), optional :: comm, fft_output
      integer :: ierr, grid(3), actual_rank, options(3), min_options(3), max_options(3)

      !! 1. Release the previous problem and select the communicator
      call ISFFT_Finalize()
      call ISFFT_Set_Communicator(comm)
      configured_planning_level = ISFFT_PLAN_MEASURE
      configured_prepare_scope = ISFFT_USE_ALL
      if (present(planning_level)) configured_planning_level = planning_level
      if (present(prepare_scope)) configured_prepare_scope = prepare_scope
      fft_output_modes = ISFFT_FFT_BOTH
      if (present(fft_output)) fft_output_modes = fft_output

      !! 2. Agree on planning scope before collective preparation
      options = [configured_planning_level, configured_prepare_scope, fft_output_modes]
      call MPI_Allreduce(options, min_options, 3, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      call MPI_Allreduce(options, max_options, 3, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      if (any(min_options /= max_options)) error stop 'planning_level and prepare_scope must agree across all MPI ranks'
      if (configured_planning_level < 0 .or. configured_planning_level > ISFFT_PLAN_MEASURE) error stop 'invalid planning_level'
      if (.not. any(configured_prepare_scope == [1, 2, 3, 4, 5, 8, 9, 16, 17, 32, 33, 63])) error stop 'invalid prepare_scope'
      if (fft_output_modes < 1 .or. fft_output_modes > 3) error stop 'invalid fft_output'
      public_fft_modes = [ibits(configured_prepare_scope, 4, 2), ibits(configured_prepare_scope, 2, 2)]
      if (configured_prepare_scope == ISFFT_USE_ALL) then
         public_fft_modes = iand(public_fft_modes, fft_output_modes)
      else if (any(iand(public_fft_modes, fft_output_modes) /= public_fft_modes)) then
         error stop 'fft_output conflicts with the selected FFT usage'
      end if
      call Policy_Configure(configured_planning_level)
      call MPI_Comm_size(ISFFT_Comm, np_size, ierr)
      call MPI_Comm_rank(ISFFT_Comm, actual_rank, ierr)
      if (dimensions < 1 .or. dimensions > 3) error stop 'invalid dimensions'
      if (any(global_grid_size < 1) .or. any(process_grid_size < 1) .or. product(process_grid_size) /= np_size .or. &
         rank /= actual_rank) error stop 'invalid grid/rank'
      if (any(process_grid_size > global_grid_size) .or. any(domain_lengths <= 0) .or. halo_width < 0 .or. &
         any(boundary_types(1:dimensions) < 0) .or. any(boundary_types(1:dimensions) > 2)) &
         error stop 'invalid geometry/boundary_types'
      if (discretization_scheme < 0 .or. discretization_scheme > 6 .or. decomposition_type < 0 .or. &
         decomposition_type > 2) &
         error stop 'invalid discretization_scheme/decomposition'
      if (dimensions == 1 .and. (np_size /= 1 .or. any(global_grid_size(2:3) /= 1) .or. decomposition_type /= 0)) &
         error stop '1D requires serial geometry'
      if (dimensions == 2 .and. &
         (global_grid_size(3) /= 1 .or. process_grid_size(3) /= 1 .or. decomposition_type == 2)) &
         error stop 'invalid 2D geometry'

      !! 3. Commit validated geometry and advance the generation
      nx_global = global_grid_size(1)
      ny_global = global_grid_size(2)
      nz_global = global_grid_size(3)
      npx0 = process_grid_size(1)
      npy0 = process_grid_size(2)
      npz0 = process_grid_size(3)
      my_id = rank
      domain_length_x = domain_lengths(1)
      domain_length_y = domain_lengths(2)
      domain_length_z = domain_lengths(3)
      configured_halo_width = halo_width
      configured_decomposition_type = decomposition_type
      configured_discretization_scheme = discretization_scheme
      configured_dimensions = dimensions
      configured_boundary_types = boundary_types
      ISFFT_plan_generation = ISFFT_plan_generation + 1
      halfspectrum_x = nx_global / 2 + 1
      halfspectrum_y = ny_global / 2 + 1
      halfspectrum_z = nz_global / 2 + 1

      !! 4. Build active topology and computation ownership
      call ISFFT_part()
      call ISFFT_Get_Internal_Process_Grid(grid(1), grid(2), grid(3))
      if (.not. btest(configured_prepare_scope, 0) .and. any(process_grid_size /= grid)) &
         error stop 'Input grid must match computation grid when redistribution is not requested'
      !! Factors are prepared on demand by the selected plans, during initialization only.

      !! 5. Prepare the requested transform and equation dependencies
      call Initialize_Plans( &
         global_grid_size, process_grid_size, grid, &
         dimensions, boundary_types, decomposition_type, &
         domain_lengths, discretization_scheme, ISFFT_plan_generation)
      if (public_fft_modes(1) /= 0) call Ensure_Public_Plan(.false.)
      if (public_fft_modes(2) /= 0) call Ensure_Public_Plan(.true.)

      !! 6. Reserve winning-plan scratch and seal execution workspaces
      call Prepare_Execution_Workspaces()
   end subroutine ! ISFFT_Initialize

   !! Call before MPI_Finalize while cached communicators and datatypes remain valid.
   !! Release borrowing plans before their transport resources and library communicator.
   subroutine ISFFT_Finalize()

      topology_ready = .false.
      if (allocated(block_cache)) deallocate (block_cache)
      if (allocated(axis_plan_cache)) deallocate (axis_plan_cache)

      call Finalize_Plans()
      call Release_Transport_Cache()

      call Release_Execution_Scratch()
      call FFT_Finalize()
      call ISFFT_Part_Finalize()
      call ISFFT_Release_Communicator()
   end subroutine ! ISFFT_Finalize

   subroutine ISFFT_Get_Local_Layout(local_size, global_start)
      integer, intent(out) :: local_size(3), global_start(3)

      if (.not. topology_ready) error stop 'ISFFT_Get_Local_Layout requires ISFFT_Initialize'
      local_size = [nx, ny, nz]
      global_start = [i_offset(npx), j_offset(npy), k_offset(npz)]
   end subroutine ! ISFFT_Get_Local_Layout

   !! Set overrides before initialization and keep them consistent across ranks.
   subroutine Set_Plan_Options(layout_choice, transport_choice)
      integer, intent(in), optional :: layout_choice, transport_choice

      if (present(layout_choice)) then
         if (layout_choice < -1 .or. layout_choice > 2) error stop 'invalid layout policy'
         layout_policy = layout_choice
      end if
      if (present(transport_choice)) then
         if (transport_choice < -1 .or. transport_choice > 4) error stop 'invalid transport policy'
         transport_policy = transport_choice
      end if
   end subroutine ! Set_Plan_Options

end module ISFFT_Initialization
