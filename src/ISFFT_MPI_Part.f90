!------------------------------------------------------------------------------
! ISFFT_MPI_Part
! Process decomposition, ownership queries and communicator topology.
!------------------------------------------------------------------------------
module ISFFT_MPI_Part
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Problem_State, only: boundary_types, decomposition_type, dimensions, i_nn, &
      i_offset, ID_XM1, ID_XP1, ID_YM1, &
      ID_YP1, ID_ZM1, ID_ZP1, ISFFT_plan_generation, &
      j_nn, j_offset, k_nn, k_offset, &
      Layout, MPI_COMM_X, MPI_COMM_Y, MPI_COMM_Z, &
      my_id, np_size, npx, npx0, &
      npy, npy0, npz, npz0, &
      nx, nx_global, ny, ny_global, &
      nz, nz_global
   implicit none
   private
   public :: Coordinates, Effective_Partition, ISFFT_Block_Decompose, ISFFT_Build_Block_Layout, ISFFT_get_id
   public :: ISFFT_Get_Internal_Process_Grid, ISFFT_mod, ISFFT_part, ISFFT_part_change, ISFFT_Part_Finalize
   public :: Axis_Position, Layout_Box, Layout_Shape, Make_Layout, Same_Layout

   type :: layout_comms
      integer :: process_grid(3), comms(3)
      type(layout_comms), pointer :: next => null()
   end type
   type(layout_comms), pointer, save :: cached_layouts => null()
   integer, save :: cache_generation = -1
   logical, save :: reuse_layout_comms = .true.

contains

   pure integer function Axis_Position(axis_order, physical_axis) result(memory_axis)
      integer, intent(in) :: axis_order(:), physical_axis
      integer :: axis

      memory_axis = 0
      do axis = 1, size(axis_order)
         if (axis_order(axis) == physical_axis) then
            memory_axis = axis
            return
         end if
      end do
   end function Axis_Position

   !------------------------------------------------------------------------------
   !! Block sizes and one-based global starts
   !------------------------------------------------------------------------------
   !! part_id is zero-based, but global_offset is a one-based global start.
   !! Layout_Box uses zero-based offsets. Convert at the interface between them.
   subroutine ISFFT_Block_Decompose(global_n, num_parts, part_id, local_n, global_offset)
      integer, intent(in) :: global_n, num_parts, part_id
      integer, intent(out) :: local_n, global_offset
      integer :: quotient, remainder

      if (global_n < 0 .or. num_parts <= 0 .or. part_id < 0 .or. part_id >= num_parts) then
         call ISFFT_Part_Abort('invalid arguments in ISFFT_Block_Decompose')
      end if

      quotient = global_n / num_parts
      remainder = mod(global_n, num_parts)
      local_n = quotient
      if (part_id < remainder) local_n = local_n + 1
      global_offset = quotient * part_id + min(part_id, remainder) + 1
   end subroutine ISFFT_Block_Decompose

   subroutine ISFFT_Build_Block_Layout(global_n, num_parts, local_sizes, global_offsets)
      integer, intent(in) :: global_n, num_parts
      integer, intent(out) :: local_sizes(0:), global_offsets(0:)
      integer :: part_id

      if (size(local_sizes) < num_parts .or. size(global_offsets) < num_parts) then
         call ISFFT_Part_Abort('layout arrays are too small in ISFFT_Build_Block_Layout')
      end if
      do part_id = 0, num_parts - 1
         call ISFFT_Block_Decompose(global_n, num_parts, part_id, local_sizes(part_id), global_offsets(part_id))
      end do
   end subroutine ISFFT_Build_Block_Layout

   !------------------------------------------------------------------------------
   !! Process-coordinate helpers
   !------------------------------------------------------------------------------
   !! This helper wraps by one period. Callers must keep i within one period
   !! of the valid process-coordinate range.
   integer function ISFFT_mod(i, n)
      use ISFFT_Problem_State

      integer, intent(in) :: i, n

      if (i .lt. 0) then
         ISFFT_mod = i + n
      else if (i .gt. n - 1) then
         ISFFT_mod = i - n
      else
         ISFFT_mod = i
      end if

   end function ISFFT_mod

   integer function ISFFT_get_id(npx1, npy1, npz1)
      use ISFFT_Problem_State

      integer :: npx1, npy1, npz1

      ISFFT_get_id = npz1 * (npx0 * npy0) + npy1 * npx0 + npx1

   end function ISFFT_get_id

   !------------------------------------------------------------------------------
   !! Active topology and communicator lifecycle
   !------------------------------------------------------------------------------
   subroutine ISFFT_part()
      use ISFFT_Problem_State

      integer :: ierr
      integer :: npx1, npy1, npz1, npx2, npy2, npz2
      integer :: active_x, active_y, active_z
      type(layout_comms), pointer :: entry
      logical :: created
      character(16) :: cache_option
      integer :: option_status

      if (np_size .ne. npx0 * npy0 * npz0) then
         call ISFFT_Part_Abort('the MPI size is not equal to npx0*npy0*npz0')
      end if
      if (max(npx0, npy0, npz0) > size(i_nn)) then
         call ISFFT_Part_Abort('the process grid exceeds the internal layout-array capacity')
      end if

      npx = mod(my_id, npx0)
      npy = mod(my_id, npx0 * npy0) / npx0
      npz = my_id / (npx0 * npy0)

      if (cache_generation /= ISFFT_plan_generation) then
         call ISFFT_Part_Finalize()
         cache_generation = ISFFT_plan_generation
         cache_option = ''
         call get_environment_variable('ISFFT_LAYOUT_CACHE', cache_option, status=option_status)
         reuse_layout_comms = .true.
         if (option_status == 0) then
            select case (trim(cache_option))
            case ('0')
               reuse_layout_comms = .false.
            case ('1', '')
               reuse_layout_comms = .true.
            case default
               call ISFFT_Part_Abort('ISFFT_LAYOUT_CACHE must be 0 or 1')
            end select
         else if (option_status /= 1) then
            call ISFFT_Part_Abort('cannot read ISFFT_LAYOUT_CACHE')
         end if
      else if (.not. reuse_layout_comms) then
         call ISFFT_Part_Finalize()
         cache_generation = ISFFT_plan_generation
      end if
      entry => cached_layouts
      do while (associated(entry))
         if (all(entry%process_grid == [npx0, npy0, npz0])) exit
         entry => entry%next
      end do
      created = .not. associated(entry)
      if (created) then
         allocate (entry)
         entry%process_grid = [npx0, npy0, npz0]
         call MPI_COMM_SPLIT(ISFFT_Comm, npz * npx0 * npy0 + npy * npx0, npx, entry%comms(1), ierr)
         call MPI_COMM_SPLIT(ISFFT_Comm, npz * npx0 * npy0 + npx, npy, entry%comms(2), ierr)
         call MPI_COMM_SPLIT(ISFFT_Comm, npy * npx0 + npx, npz, entry%comms(3), ierr)
         entry%next => cached_layouts
         cached_layouts => entry
      end if
      MPI_COMM_X = entry%comms(1)
      MPI_COMM_Y = entry%comms(2)
      MPI_COMM_Z = entry%comms(3)

      call ISFFT_Build_Block_Layout(nx_global, npx0, i_nn, i_offset)
      call ISFFT_Build_Block_Layout(ny_global, npy0, j_nn, j_offset)
      call ISFFT_Build_Block_Layout(nz_global, npz0, k_nn, k_offset)
      call ISFFT_Validate_Block_Layout(nx_global, npx0, i_nn, i_offset, 'x')
      call ISFFT_Validate_Block_Layout(ny_global, npy0, j_nn, j_offset, 'y')
      call ISFFT_Validate_Block_Layout(nz_global, npz0, k_nn, k_offset, 'z')
      nx = i_nn(npx)
      ny = j_nn(npy)
      nz = k_nn(npz)

      active_x = min(nx_global, npx0)
      if (npx < active_x) then
         npx1 = ISFFT_mod(npx - 1, active_x)
         npx2 = ISFFT_mod(npx + 1, active_x)
         ID_XM1 = npz * (npx0 * npy0) + npy * npx0 + npx1
         ID_XP1 = npz * (npx0 * npy0) + npy * npx0 + npx2
         if (boundary_types(1) .ne. 0 .and. npx .eq. 0) ID_XM1 = MPI_PROC_NULL
         if (boundary_types(1) .ne. 0 .and. npx .eq. active_x - 1) ID_XP1 = MPI_PROC_NULL
      else
         ID_XM1 = MPI_PROC_NULL
         ID_XP1 = MPI_PROC_NULL
      end if

      active_y = min(ny_global, npy0)
      if (npy < active_y) then
         npy1 = ISFFT_mod(npy - 1, active_y)
         npy2 = ISFFT_mod(npy + 1, active_y)
         ID_YM1 = npz * (npx0 * npy0) + npy1 * npx0 + npx
         ID_YP1 = npz * (npx0 * npy0) + npy2 * npx0 + npx
         if (boundary_types(2) .ne. 0 .and. npy .eq. 0) ID_YM1 = MPI_PROC_NULL
         if (boundary_types(2) .ne. 0 .and. npy .eq. active_y - 1) ID_YP1 = MPI_PROC_NULL
      else
         ID_YM1 = MPI_PROC_NULL
         ID_YP1 = MPI_PROC_NULL
      end if

      active_z = min(nz_global, npz0)
      if (npz < active_z) then
         npz1 = ISFFT_mod(npz - 1, active_z)
         npz2 = ISFFT_mod(npz + 1, active_z)
         ID_ZM1 = npz1 * (npx0 * npy0) + npy * npx0 + npx
         ID_ZP1 = npz2 * (npx0 * npy0) + npy * npx0 + npx
         if (boundary_types(3) .ne. 0 .and. npz .eq. 0) ID_ZM1 = MPI_PROC_NULL
         if (boundary_types(3) .ne. 0 .and. npz .eq. active_z - 1) ID_ZP1 = MPI_PROC_NULL
      else
         ID_ZM1 = MPI_PROC_NULL
         ID_ZP1 = MPI_PROC_NULL
      end if

      !! A cached layout switch updates local metadata only; no synchronization
      !! is required. Keep the original barrier on first construction.
      if (created) call MPI_barrier(ISFFT_Comm, ierr)

   end subroutine ISFFT_part

   !! This updates active topology metadata. Field movement is handled by callers.
   subroutine ISFFT_part_change(npx0_new, npy0_new, npz0_new)
      use ISFFT_Problem_State
      integer :: npx0_new, npy0_new, npz0_new

      npx0 = npx0_new

      npy0 = npy0_new

      npz0 = npz0_new

      if (np_size .ne. npx0 * npy0 * npz0) then
         call ISFFT_Part_Abort('the changed process grid does not match the MPI size')
      end if

      call ISFFT_part

   end subroutine ISFFT_part_change

   subroutine ISFFT_Part_Finalize()
      type(layout_comms), pointer :: entry
      integer :: axis, ierr

      do while (associated(cached_layouts))
         entry => cached_layouts
         cached_layouts => entry%next
         do axis = 1, 3
            call MPI_Comm_free(entry%comms(axis), ierr)
         end do
         deallocate (entry)
      end do
      MPI_COMM_X = MPI_COMM_NULL
      MPI_COMM_Y = MPI_COMM_NULL
      MPI_COMM_Z = MPI_COMM_NULL
      cache_generation = -1
   end subroutine ISFFT_Part_Finalize

   !------------------------------------------------------------------------------
   !! Ownership validation and collective abort
   !------------------------------------------------------------------------------
   subroutine ISFFT_Validate_Block_Layout(global_n, num_parts, local_sizes, global_offsets, axis_name)
      integer, intent(in) :: global_n, num_parts
      integer, intent(in) :: local_sizes(0:), global_offsets(0:)
      character(len=*), intent(in) :: axis_name
      integer :: part_id, expected_size, expected_offset, total_size

      total_size = 0
      do part_id = 0, num_parts - 1
         call ISFFT_Block_Decompose(global_n, num_parts, part_id, expected_size, expected_offset)
         if (local_sizes(part_id) /= expected_size .or. global_offsets(part_id) /= expected_offset) then
            call ISFFT_Part_Abort('invalid canonical block layout on axis '//trim(axis_name))
         end if
         total_size = total_size + local_sizes(part_id)
      end do
      if (total_size /= global_n) then
         call ISFFT_Part_Abort('incomplete canonical block layout on axis '//trim(axis_name))
      end if
   end subroutine ISFFT_Validate_Block_Layout

   subroutine ISFFT_Part_Abort(message)
      character(len=*), intent(in) :: message
      integer :: ierr

      if (my_id == 0) write (*, '(A,A)') 'ERROR: ', trim(message)
      call MPI_Abort(ISFFT_Comm, 1, ierr)
      stop 1
   end subroutine ISFFT_Part_Abort

   !------------------------------------------------------------------------------
   !! Layout geometry and zero-based box offsets
   !------------------------------------------------------------------------------
   pure function Effective_Partition(field_layout) result(domain)
      type(Layout), intent(in) :: field_layout
      integer :: domain(3)

      !! Layout_Box interprets nonpositive partition entries as global lengths.
      domain = merge(field_layout%partition_shape, field_layout%global_shape, field_layout%partition_shape > 0)
   end function ! Effective_Partition

   logical function Same_Layout(left_layout, right_layout)
      type(Layout), intent(in) :: left_layout, right_layout

      Same_Layout = all(left_layout%global_shape == right_layout%global_shape) .and. &
         all(left_layout%process_grid == right_layout%process_grid) .and. &
         all(left_layout%axis_order == right_layout%axis_order) .and. &
         all(Effective_Partition(left_layout) == Effective_Partition(right_layout)) .and. &
         left_layout%replicated_axis == right_layout%replicated_axis .and. &
         (left_layout%complex_data .eqv. right_layout%complex_data)
      if (allocated(left_layout%rank_boxes) .neqv. allocated(right_layout%rank_boxes)) then
         Same_Layout = .false.
         return
      end if
      if (allocated(left_layout%rank_boxes) .and. allocated(right_layout%rank_boxes)) then
         if (any(shape(left_layout%rank_boxes) /= shape(right_layout%rank_boxes))) then
            Same_Layout = .false.
            return
         end if
         Same_Layout = Same_Layout .and. all(left_layout%rank_boxes == right_layout%rank_boxes)
      end if
   end function ! Same_Layout

   function Make_Layout(global_shape, process_grid, axis_order, complex_data, replicated_axis) result(l)
      integer, intent(in) :: global_shape(3), process_grid(3), axis_order(3)
      logical, intent(in) :: complex_data
      integer, intent(in), optional :: replicated_axis
      type(Layout) :: l

      l%global_shape = global_shape
      l%process_grid = process_grid
      l%axis_order = axis_order
      l%complex_data = complex_data
      if (present(replicated_axis)) l%replicated_axis = replicated_axis
   end function ! Make_Layout

   subroutine Coordinates(field_layout, rank, process_coordinates)
      type(Layout), intent(in) :: field_layout
      integer, intent(in) :: rank
      integer, intent(out) :: process_coordinates(3)

      process_coordinates = &
         [mod(rank, field_layout%process_grid(1)), &
         mod(rank / field_layout%process_grid(1), field_layout%process_grid(2)), rank / &
         (field_layout%process_grid(1) * field_layout%process_grid(2))]
   end subroutine ! Coordinates

   subroutine Layout_Box(field_layout, rank, local_shape, global_offset)
      type(Layout), intent(in) :: field_layout
      integer, intent(in) :: rank
      integer, intent(out) :: local_shape(3), global_offset(3)
      integer :: c(3), a, q, r, domain

      if (allocated(field_layout%rank_boxes)) then
         local_shape = field_layout%rank_boxes(:3, rank)
         global_offset = field_layout%rank_boxes(4:, rank)
         return
      end if
      call Coordinates(field_layout, rank, c)
      do a = 1, 3
         domain = merge(field_layout%partition_shape(a), field_layout%global_shape(a), &
            field_layout%partition_shape(a) > 0)
         q = domain / field_layout%process_grid(a)
         r = mod(domain, field_layout%process_grid(a))
         local_shape(a) = q
         if (c(a) < r) local_shape(a) = local_shape(a) + 1
         global_offset(a) = q * c(a) + min(c(a), r)
         local_shape(a) = max(0, min(local_shape(a), field_layout%global_shape(a) - global_offset(a)))
      end do
      if (field_layout%replicated_axis > 0) then
         local_shape(field_layout%replicated_axis) = field_layout%global_shape(field_layout%replicated_axis)
         global_offset(field_layout%replicated_axis) = 0
      end if
   end subroutine ! Layout_Box

   function Layout_Shape(field_layout, rank) result(n)
      type(Layout), intent(in) :: field_layout
      integer, intent(in) :: rank
      integer :: n(3), off(3), physical(3)

      call Layout_Box(field_layout, rank, physical, off)
      n = physical(field_layout%axis_order)
   end function ! Layout_Shape

   !------------------------------------------------------------------------------
   !! Computation process-grid selection
   !------------------------------------------------------------------------------
   subroutine ISFFT_Get_Internal_Process_Grid(processes_x, processes_y, processes_z)
      integer, intent(out) :: processes_x, processes_y, processes_z

      select case (dimensions)
      case (1)
         processes_x = 1
         processes_y = 1
         processes_z = 1
      case (2)
         select case (decomposition_type)
         case (0)
            processes_x = npx0
            processes_y = npy0
            processes_z = 1
         case (1)
            processes_x = 1
            processes_y = np_size
            processes_z = 1
         case default
            call Redistribution_Abort('dimensions=2 does not support pencil decomposition')
         end select
      case (3)
         select case (decomposition_type)
         case (0)
            processes_x = npx0
            processes_y = npy0
            processes_z = npz0
         case (1)
            processes_x = 1
            processes_y = 1
            processes_z = np_size
            if (np_size > min(nx_global, ny_global, nz_global)) then
               call Redistribution_Abort('The 3D slab decomposition is not suitable')
            end if
         case (2)
            call Select_Pencil_Process_Grid(np_size, processes_y, processes_z)
            processes_x = 1
         case default
            call Redistribution_Abort('dimensions=3 requires block, slab, or pencil decomposition')
         end select
      case default
         call Redistribution_Abort('redistribution requires dimensions=1, 2, or 3')
      end select

   end subroutine ISFFT_Get_Internal_Process_Grid

   subroutine Select_Pencil_Process_Grid(num_processes, num_processes_y, num_processes_z)
      integer, intent(in) :: num_processes
      integer, intent(out) :: num_processes_y, num_processes_z
      integer :: divisor

      num_processes_y = num_processes
      num_processes_z = 1
      do divisor = int(sqrt(real(num_processes))), 1, -1
         if (mod(num_processes, divisor) == 0) then
            num_processes_z = divisor
            num_processes_y = num_processes / divisor
            return
         end if
      end do
   end subroutine Select_Pencil_Process_Grid

   !------------------------------------------------------------------------------
   !! Redistribution preconditions and diagnostics
   !------------------------------------------------------------------------------
   subroutine Redistribution_Abort(message)
      character(len=*), intent(in) :: message
      integer :: ierr

      if (my_id == 0) write (*, '(A,A)') 'ERROR: ', trim(message)
      call MPI_Abort(ISFFT_Comm, 1, ierr)
      stop 1
   end subroutine Redistribution_Abort

end module ISFFT_MPI_Part
