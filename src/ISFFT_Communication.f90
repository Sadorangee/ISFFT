!------------------------------------------------------------------------------
! ISFFT_Communication
! Transpose and redistribution mappings, resources and data movement.
!------------------------------------------------------------------------------
module ISFFT_Communication
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_COMPLEX_MPI_TYPE, ISFFT_INT64_KIND, ISFFT_REAL_KIND, ISFFT_REAL_MPI_TYPE
   use ISFFT_Plans, only: Axis_Plan, default_chunk_elements, default_stream_elements, Kernel_Plan, &
      redistribution_forward_plan, redistribution_inverse_plan, Transfer_Plan, Transform_Plan
   use ISFFT_Problem_State, only: halo_width, internal_layout, Layout, my_id, &
      npx0, npy0, npz0, nx, &
      ny, nz, original_layout
   use ISFFT_Workspace, only: block_counts, block_displs, block_recv_r, block_recv_z, &
      block_send_r, block_send_z, block_tile_r, block_tile_z, &
      Bound_Layout_Workspace, Get_Complex, Get_Real, Park_Field, &
      recv_complex_buffer, recv_real_buffer, redistribution_field, Reserve_R, &
      Reserve_Z, send_complex_buffer, send_real_buffer, Spectral_Field
   use ISFFT_MPI_Part, only: Coordinates, Effective_Partition, ISFFT_part_change, Layout_Box, &
      Layout_Shape, Same_Layout
   implicit none
   private
   public :: Begin_Block_Transfer, Block_Transfer, Build_Block_Gather, Build_Transfer, datatype_a2a
   public :: datatype_p2p, End_Block_Transfer, Gather_Block_Tile, memorysave_p2p, Move_Field
   public :: Move_State, packed_a2a, packed_p2p, Redistribute_Field, Release_Transfer
   public :: Release_Transform_Plan, Release_Transport_Cache, Require_Redistribution_Initialization, &
      resources, Select_Transfer_Method
   public :: Set_Explicit_Boxes, Store_Block_Tile, Transfer_Method_Name, Transport_Resource, Transport_Resource_Count

   integer, parameter :: packing_tile_width = 16
   integer, parameter :: packed_a2a = 0, datatype_a2a = 1, packed_p2p = 2, datatype_p2p = 3, memorysave_p2p = 4

   type :: Transport_Resource
      type(Transfer_Plan) :: geometry
      integer :: comm = MPI_COMM_NULL, transfer_kind = 0
      logical :: packed_p2p_active(0:1) = .false.
      logical :: tuned(0:1) = .false., direct(0:1) = .false., streamed(0:1) = .false., sparse(0:1) = .false.
      integer(ISFFT_INT64_KIND) :: width(0:1) = 0, limit = 0
   end type
   type(Transport_Resource), allocatable, save :: resources(:)

   interface Move_Field
      module procedure Move_Real, Move_Complex
   end interface

   type :: Block_Transfer
      type(Layout) :: source, destination
      type(Spectral_Field) :: result
      integer :: comm, parts, axis, nin(3), nout(3), oi(3), oo(3), others(2)
      integer :: lines, capacity, local_n, full_n, out_n
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: tile_r(:, :, :) => null()
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: tile_z(:, :, :) => null()
      real(kind=ISFFT_REAL_KIND), pointer :: sr(:) => null(), rr(:) => null()
      complex(kind=ISFFT_REAL_KIND), pointer :: sz(:) => null(), rz(:) => null()
   end type Block_Transfer

contains

   !------------------------------------------------------------------------------
   !! Transport resource registration and reuse
   !------------------------------------------------------------------------------
   integer function Transport_Resource_Count() result(n)

      n = 0
      if (allocated(resources)) n = size(resources)
   end function ! Transport_Resource_Count

   subroutine Build_Transfer(plan, source_layout, destination_layout)
      type(Transfer_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: source_layout, destination_layout

      call Prepare_Transport(plan, source_layout, destination_layout, 0)
   end subroutine ! Build_Transfer

   subroutine Build_Block_Gather(plan, source_layout, destination_layout)
      type(Transfer_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: source_layout, destination_layout

      call Prepare_Transport(plan, source_layout, destination_layout, 1)
   end subroutine ! Build_Block_Gather

   subroutine Prepare_Transport(plan, source_layout, destination_layout, transfer_kind)
      type(Transfer_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: source_layout, destination_layout
      type(Layout) :: normalized_source, normalized_destination
      integer, intent(in) :: transfer_kind
      type(Transport_Resource), allocatable :: grown_resources(:)
      integer :: resource_index, cached_index, shared_cached_index, ierr, local_rank, status, resource_count, rank_count
      integer(ISFFT_INT64_KIND) :: chunk_capacity
      character(64) :: chunk_option

      normalized_source = source_layout
      normalized_destination = destination_layout
      normalized_source%partition_shape = Effective_Partition(source_layout)
      normalized_destination%partition_shape = Effective_Partition(destination_layout)
      call Release_Transfer(plan)
      if (transfer_kind == 0 .and. Same_Layout(normalized_source, normalized_destination)) then
         call MPI_Comm_size(ISFFT_Comm, rank_count, ierr)
         if (product(normalized_source%process_grid) /= rank_count) error stop 'invalid identity layout grid'
         plan%source = normalized_source
         plan%destination = normalized_destination
         plan%local_only = .true.
         plan%ready = .true.
         return
      end if
      !! A single rank needs local layout metadata, never a transport search.
      !! Keep the same local execution path, including nontrivial axis permutations.
      call MPI_Comm_size(ISFFT_Comm, rank_count, ierr)
      if (rank_count == 1) then
         call Build_Transfer_Raw(plan, normalized_source, normalized_destination)
         return
      end if
      chunk_capacity = default_chunk_elements
      call MPI_Comm_rank(ISFFT_Comm, local_rank, ierr)
      if (local_rank == 0) then
         call get_environment_variable('ISFFT_MPI_CHUNK_ELEMENTS', chunk_option, status=status)
         if (status == 0 .and. len_trim(chunk_option) > 0) then
            read (chunk_option, *, iostat=status) chunk_capacity
            if (status /= 0 .or. chunk_capacity < 1 .or. chunk_capacity > int(huge(0), ISFFT_INT64_KIND)) &
               error stop 'Invalid chunk capacity'
         end if
      end if
      call MPI_Bcast(chunk_capacity, 1, MPI_INTEGER8, 0, ISFFT_Comm, ierr)
      cached_index = 0
      if (allocated(resources)) then
         do resource_index = 1, size(resources)
            if (resources(resource_index)%comm /= ISFFT_Comm .or. &
               resources(resource_index)%transfer_kind /= transfer_kind .or. &
               resources(resource_index)%limit /= chunk_capacity) cycle
            if (Same_Layout(normalized_source, resources(resource_index)%geometry%source) .and. &
               Same_Layout(normalized_destination, resources(resource_index)%geometry%destination)) cached_index = &
               resource_index
         end do
      end if
      call MPI_Allreduce(cached_index, shared_cached_index, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      if (shared_cached_index > 0) then
         plan = resources(cached_index)%geometry
         plan%resource_id = cached_index
         return
      end if
      if (transfer_kind == 0) then
         call Build_Transfer_Raw(plan, normalized_source, normalized_destination)
      else
         call Build_Block_Gather_Raw(plan, normalized_source, normalized_destination)
      end if
      if (.not. allocated(resources)) then
         allocate (resources(0))
      end if
      resource_count = size(resources)
      allocate (grown_resources(resource_count + 1))
      grown_resources(:resource_count) = resources
      grown_resources(resource_count + 1)%geometry = plan
      grown_resources(resource_count + 1)%comm = ISFFT_Comm
      grown_resources(resource_count + 1)%transfer_kind = transfer_kind
      grown_resources(resource_count + 1)%limit = chunk_capacity
      call move_alloc(grown_resources, resources)
      plan%resource_id = resource_count + 1
   end subroutine ! Prepare_Transport

   !! Release borrowing plans before the cache. Cached geometries own MPI handles.
   subroutine Release_Transport_Cache()
      integer :: i

      if (.not. allocated(resources)) return
      do i = 1, size(resources)
         call Release_Transfer(resources(i)%geometry)
      end do
      deallocate (resources)
   end subroutine ! Release_Transport_Cache

   !! Explicit backend-owned boxes use the same transport algorithms and cache.
   !! They describe communication geometry only, not a new FFT decomposition.
   subroutine Set_Explicit_Boxes(field_layout, local_shape, global_offset)
      type(Layout), intent(inout) :: field_layout
      integer, intent(in) :: local_shape(3), global_offset(3)
      integer :: np, ierr, i, j, lo(3), hi(3)
      integer(ISFFT_INT64_KIND) :: volume

      call MPI_Comm_size(ISFFT_Comm, np, ierr)
      if (allocated(field_layout%rank_boxes)) deallocate (field_layout%rank_boxes)
      allocate (field_layout%rank_boxes(6, 0:np - 1))
      call MPI_Allgather([local_shape, global_offset], 6, MPI_INTEGER, field_layout%rank_boxes, 6, MPI_INTEGER, &
         ISFFT_Comm, ierr)
      volume = 0
      do i = 0, np - 1
         if (any(field_layout%rank_boxes(:, i) < 0) .or. &
            any(field_layout%rank_boxes(:3, i) + field_layout%rank_boxes(4:, i) > field_layout%global_shape)) &
            error stop 'Invalid explicit box bounds'
         volume = volume + product(int(field_layout%rank_boxes(:3, i), ISFFT_INT64_KIND))
         do j = 0, i - 1
            lo = max(field_layout%rank_boxes(4:, i), field_layout%rank_boxes(4:, j))
            hi = min(field_layout%rank_boxes(4:, i) + field_layout%rank_boxes(:3, i), &
               field_layout%rank_boxes(4:, j) + field_layout%rank_boxes(:3, j))
            if (all(hi > lo)) error stop 'Overlapping explicit ownership boxes'
         end do
      end do
      if (volume /= product(int(field_layout%global_shape, ISFFT_INT64_KIND))) &
         error stop 'Explicit boxes do not cover global domain'
   end subroutine ! Set_Explicit_Boxes

   subroutine Intersection(source_layout, destination_layout, source_rank, destination_rank, intersection_start, &
      intersection_shape)
      type(Layout), intent(in) :: source_layout, destination_layout
      integer, intent(in) :: source_rank, destination_rank
      integer, intent(out) :: intersection_start(3), intersection_shape(3)
      integer :: ns(3), nd(3), os(3), od(3), cs(3), cd(3), a

      call Layout_Box(source_layout, source_rank, ns, os)
      call Layout_Box(destination_layout, destination_rank, nd, od)
      intersection_start = max(os, od)
      intersection_shape = max(0, min(os + ns, od + nd) - intersection_start)
      !! A replicated source contributes exactly one copy. Preserve the receiver's
      !! replica coordinate when possible, so contraction is a purely local slice.
      a = source_layout%replicated_axis
      if (a > 0) then
         call Coordinates(source_layout, source_rank, cs)
         call Coordinates(destination_layout, destination_rank, cd)
         if (cs(a) /= mod(cd(a), source_layout%process_grid(a))) intersection_shape = 0
      end if
   end subroutine ! Intersection

   !------------------------------------------------------------------------------
   !! Raw transport construction and MPI handle ownership
   !------------------------------------------------------------------------------
   subroutine Build_Block_Gather_Raw(plan, source_layout, destination_layout)
      type(Transfer_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: source_layout, destination_layout
      integer :: rank, np, ierr, axis, c(3), peer(3), color

      call Release_Transfer(plan)
      axis = destination_layout%replicated_axis
      if (axis < 1 .or. axis > 3) error stop 'block gather requires a replicated target axis'
      if (source_layout%replicated_axis /= 0 .or. &
         any(source_layout%process_grid /= destination_layout%process_grid) .or. &
         any(Effective_Partition(source_layout) /= Effective_Partition(destination_layout)) .or. &
         any(source_layout%global_shape /= destination_layout%global_shape) .or. &
         (source_layout%complex_data .neqv. destination_layout%complex_data)) &
         error stop 'incompatible block gather layouts'
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call MPI_Comm_size(ISFFT_Comm, np, ierr)
      if (product(source_layout%process_grid) /= np) error stop 'invalid block gather grid'
      plan%source = source_layout
      plan%destination = destination_layout
      plan%local_only = .false.
      call Coordinates(source_layout, rank, c)
      peer = c
      peer(axis) = 0
      color = peer(1) + source_layout%process_grid(1) * (peer(2) + source_layout%process_grid(2) * peer(3))
      call MPI_Comm_split(ISFFT_Comm, color, c(axis), plan%gather_comm, ierr)
   !! Execute_Block_Axis computes bounded tile counts itself. This owns only
   !! its gather topology, not a ready whole-field Move_Field plan.
   end subroutine ! Build_Block_Gather_Raw

   subroutine Build_Transfer_Raw(plan, source_layout, destination_layout)
      type(Transfer_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: source_layout, destination_layout
      integer :: rank, np, ierr, q, ns(3), nd(3), os(3), od(3), lo(3), base, bytes, local_flag, global_flag
      integer :: axis, c(3), peer(3), color, j, peer_rank, group_size, color_stride, status
      integer(ISFFT_INT64_KIND) :: local_max
      character(64) :: option

      call Release_Transfer(plan)
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call MPI_Comm_size(ISFFT_Comm, np, ierr)
      if (any(source_layout%global_shape /= destination_layout%global_shape) .or. &
         (source_layout%complex_data .neqv. destination_layout%complex_data)) &
         error stop 'transfer requires equal global shapes and data type'
      if (product(source_layout%process_grid) /= np .or. product(destination_layout%process_grid) /= np) &
         error stop 'invalid layout grid'
      plan%source = source_layout
      plan%destination = destination_layout
      !! Identity is an execution no-op (or local copy through Move_Field).
      !! No peer counts, derived datatypes or communication candidates are needed.
      if (Same_Layout(source_layout, destination_layout)) then
         plan%local_only = .true.
         plan%chunked = .false.
         plan%max_peer = 0
         plan%ready = .true.
         return
      end if
      allocate (plan%send_counts(0:np - 1), plan%recv_counts(0:np - 1), plan%send_displs(0:np - 1), &
         plan%recv_displs(0:np - 1))
      allocate (plan%send_box_offsets(3, 0:np - 1), plan%recv_box_offsets(3, 0:np - 1), &
         plan%send_box_sizes(3, 0:np - 1), &
         plan%recv_box_sizes(3, 0:np - 1))
      allocate (plan%send_types(0:np - 1), plan%recv_types(0:np - 1), plan%send_type_counts(0:np - 1), &
         plan%recv_type_counts(0:np - 1), plan%send_byte_displs(0:np - 1), plan%recv_byte_displs(0:np - 1))
      base = ISFFT_REAL_MPI_TYPE
      if (source_layout%complex_data) call MPI_Type_match_size( &
         MPI_TYPECLASS_COMPLEX, 2 * storage_size(0.0_ISFFT_REAL_KIND) / 8, base, ierr)
      call MPI_Type_size(base, bytes, ierr)
      plan%send_types = base
      plan%recv_types = base
      plan%send_type_counts = 0
      plan%recv_type_counts = 0
      plan%send_byte_displs = 0
      plan%recv_byte_displs = 0
      call Layout_Box(source_layout, rank, ns, os)
      call Layout_Box(destination_layout, rank, nd, od)
      do q = 0, np - 1
         call Intersection(source_layout, destination_layout, rank, q, lo, plan%send_box_sizes(:, q))
         plan%send_box_offsets(:, q) = lo - os
         call Intersection(source_layout, destination_layout, q, rank, lo, plan%recv_box_sizes(:, q))
         plan%recv_box_offsets(:, q) = lo - od
         plan%send_counts(q) = product(int(plan%send_box_sizes(:, q), ISFFT_INT64_KIND))
         plan%recv_counts(q) = product(int(plan%recv_box_sizes(:, q), ISFFT_INT64_KIND))
         plan%send_displs(q) = 0_ISFFT_INT64_KIND
         plan%recv_displs(q) = 0_ISFFT_INT64_KIND
         if (q > 0) then
            plan%send_displs(q) = plan%send_displs(q - 1) + plan%send_counts(q - 1)
            plan%recv_displs(q) = plan%recv_displs(q - 1) + plan%recv_counts(q - 1)
         end if
      end do
      local_flag = 1
      do q = 0, np - 1
         if (q /= rank .and. (plan%send_counts(q) > 0 .or. plan%recv_counts(q) > 0)) local_flag = 0
      end do
      global_flag = local_flag
      if (np > 1) call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      plan%local_only = global_flag == 1
      if (plan%local_only) then
         !! Local axis-order changes use only the stride/offset metadata above.
         !! They must not prepare MPI datatypes or a communication trial.
         plan%ready = .true.
         return
      end if
      do q = 0, np - 1
         if (plan%send_counts(q) > 0) then
            call Box_Type( &
               ns, source_layout%axis_order, plan%send_box_offsets(:, q), plan%send_box_sizes(:, q), &
               base, bytes, &
               plan%send_types(q), plan%send_byte_displs(q))
            plan%send_type_counts(q) = 1
         end if
         if (plan%recv_counts(q) > 0) then
            call Box_Type( &
               nd, destination_layout%axis_order, plan%recv_box_offsets(:, q), plan%recv_box_sizes(:, q), &
               base, bytes, &
               plan%recv_types(q), plan%recv_byte_displs(q))
            plan%recv_type_counts(q) = 1
         end if
      end do
      axis = destination_layout%replicated_axis
      if (axis > 0 .and. source_layout%replicated_axis == 0 .and. &
         all(source_layout%process_grid == destination_layout%process_grid) .and. &
         all(Effective_Partition(source_layout) == Effective_Partition(destination_layout))) then
         call Coordinates(source_layout, rank, c)
         peer = c
         peer(axis) = 0
         color = peer(1) + source_layout%process_grid(1) * (peer(2) + source_layout%process_grid(2) * peer(3))
         call MPI_Comm_split(ISFFT_Comm, color, c(axis), plan%gather_comm, ierr)
         allocate (plan%gather_counts(0:source_layout%process_grid(axis) - 1), &
            plan%gather_displs(0:source_layout%process_grid(axis) - 1))
         do j = 0, source_layout%process_grid(axis) - 1
            peer = c
            peer(axis) = j
            peer_rank = peer(1) + source_layout%process_grid(1) * (peer(2) + source_layout%process_grid(2) * peer(3))
            plan%gather_counts(j) = plan%recv_counts(peer_rank)
            plan%gather_displs(j) = plan%recv_displs(peer_rank)
         end do
         plan%gather_count = product(int(ns, ISFFT_INT64_KIND))
      end if
      if (.not. plan%local_only .and. plan%gather_comm == MPI_COMM_NULL) then
         call Coordinates(source_layout, rank, c)
         color = 0
         color_stride = 1
         do axis = 1, 3
            if (allocated(source_layout%rank_boxes) .or. allocated(destination_layout%rank_boxes)) cycle
            if (source_layout%process_grid(axis) /= destination_layout%process_grid(axis) .or. &
               axis == source_layout%replicated_axis .or. axis == destination_layout%replicated_axis) cycle
            if (product(source_layout%process_grid(:axis - 1)) /= product(destination_layout%process_grid(:axis - &
               1))) cycle
            if (merge(source_layout%partition_shape(axis), source_layout%global_shape(axis), &
               source_layout%partition_shape(axis) > 0) /= &
               merge(destination_layout%partition_shape(axis), destination_layout%global_shape(axis), &
               destination_layout%partition_shape(axis) > 0)) cycle
            color = color + color_stride * c(axis)
            color_stride = color_stride * source_layout%process_grid(axis)
         end do
         call MPI_Comm_split(ISFFT_Comm, color, rank, plan%exchange_comm, ierr)
         call MPI_Comm_size(plan%exchange_comm, group_size, ierr)
         allocate (plan%peers(group_size))
         call MPI_Allgather(rank, 1, MPI_INTEGER, plan%peers, 1, MPI_INTEGER, plan%exchange_comm, ierr)
         allocate (plan%sparse_send(count(plan%send_counts(plan%peers) > 0)), &
            plan%sparse_recv(count(plan%recv_counts(plan%peers) > 0)))
         plan%sparse_send = pack([(j - 1, j=1, group_size)], plan%send_counts(plan%peers) > 0)
         plan%sparse_recv = pack([(j - 1, j=1, group_size)], plan%recv_counts(plan%peers) > 0)
         allocate (plan%requests(size(plan%sparse_send) + size(plan%sparse_recv)))
         plan%requests = MPI_REQUEST_NULL
      end if
      !! Decide the common communication path once, outside execution timing.
      plan%chunk_limit = default_chunk_elements
      if (rank == 0) then
         call get_environment_variable('ISFFT_MPI_CHUNK_ELEMENTS', option, status=status)
         if (status == 0 .and. len_trim(option) > 0) then
            read (option, *, iostat=status) plan%chunk_limit
            if (status /= 0 .or. plan%chunk_limit < 1 .or. &
               plan%chunk_limit > int(huge(0), ISFFT_INT64_KIND)) error stop 'Invalid ISFFT_MPI_CHUNK_ELEMENTS'
         end if
      end if
      call MPI_Bcast(plan%chunk_limit, 1, MPI_INTEGER8, 0, ISFFT_Comm, ierr)
      local_max = max(maxval(plan%send_counts), maxval(plan%recv_counts), plan%gather_count)
      call MPI_Allreduce(local_max, plan%max_peer, 1, MPI_INTEGER8, MPI_MAX, ISFFT_Comm, ierr)
      local_flag = 0
      if (local_max > plan%chunk_limit .or. &
         maxval(plan%send_displs) > int(huge(0), ISFFT_INT64_KIND) .or. maxval(plan%recv_displs) &
         > int(huge(0), ISFFT_INT64_KIND)) local_flag = 1
      call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      plan%chunked = global_flag == 1
      plan%streamed = plan%chunked
      plan%stream_width = min(plan%chunk_limit, default_stream_elements)
      plan%ready = .true.
   end subroutine ! Build_Transfer_Raw

   subroutine Box_Type( &
      physical_shape, axis_order, box_offset, box_size, &
      base_datatype, bytes_per_element, &
      datatype, byte_displacement)
      integer, intent(in) :: physical_shape(3), axis_order(3), box_offset(3), box_size(3), base_datatype, &
         bytes_per_element
      integer, intent(out) :: datatype, byte_displacement
      integer :: line_datatype, plane_datatype, box_datatype, ierr
      integer(MPI_ADDRESS_KIND) :: element_strides(3), byte_offset(1)

      element_strides(axis_order(1)) = 1
      element_strides(axis_order(2)) = int(physical_shape(axis_order(1)), MPI_ADDRESS_KIND)
      element_strides(axis_order(3)) = int(physical_shape(axis_order(1)), MPI_ADDRESS_KIND) * &
         int(physical_shape(axis_order(2)), MPI_ADDRESS_KIND)
      byte_offset(1) = sum(int(box_offset, MPI_ADDRESS_KIND) * element_strides) * &
         int(bytes_per_element, MPI_ADDRESS_KIND)
      !! MPI-3 Alltoallw has int byte displacements. Embed the full address-sized offset in the datatype.
      byte_displacement = 0
      call MPI_Type_create_hvector(box_size(1), 1, element_strides(1) * bytes_per_element, base_datatype, &
         line_datatype, ierr)
      call MPI_Type_create_hvector(box_size(2), 1, element_strides(2) * bytes_per_element, line_datatype, &
         plane_datatype, ierr)
      call MPI_Type_create_hvector(box_size(3), 1, element_strides(3) * bytes_per_element, plane_datatype, &
         box_datatype, ierr)
      call MPI_Type_create_struct(1, [1], byte_offset, [box_datatype], datatype, ierr)
      call MPI_Type_commit(datatype, ierr)
      call MPI_Type_free(line_datatype, ierr)
      call MPI_Type_free(plane_datatype, ierr)
      call MPI_Type_free(box_datatype, ierr)
   end subroutine ! Box_Type

   !! resource_id > 0 marks a borrower. The cache releases the shared MPI handles.
   subroutine Release_Transfer(plan)
      type(Transfer_Plan), intent(inout) :: plan
      type(Transfer_Plan) :: empty_plan
      integer :: q, ierr

      if (plan%resource_id > 0) then
         plan = empty_plan
         return
      end if
      if (allocated(plan%send_types)) then
         do q = 0, size(plan%send_types) - 1
            if (plan%send_type_counts(q) > 0) call MPI_Type_free(plan%send_types(q), ierr)
            if (plan%recv_type_counts(q) > 0) call MPI_Type_free(plan%recv_types(q), ierr)
         end do
         deallocate (plan%send_counts, plan%recv_counts, plan%send_displs, plan%recv_displs, plan%send_box_offsets, &
            plan%recv_box_offsets, plan%send_box_sizes, plan%recv_box_sizes, plan%send_types, plan%recv_types, &
            plan%send_type_counts, &
            plan%recv_type_counts, plan%send_byte_displs, plan%recv_byte_displs)
      end if
      if (plan%exchange_comm /= MPI_COMM_NULL) call MPI_Comm_free(plan%exchange_comm, ierr)
      plan%exchange_comm = MPI_COMM_NULL
      if (allocated(plan%peers)) deallocate (plan%peers)
      if (allocated(plan%sparse_send)) deallocate (plan%sparse_send, plan%sparse_recv, plan%requests)
      if (plan%gather_comm /= MPI_COMM_NULL) call MPI_Comm_free(plan%gather_comm, ierr)
      if (allocated(plan%gather_counts)) deallocate (plan%gather_counts, plan%gather_displs)
      plan%gather_comm = MPI_COMM_NULL
      plan%gather_count = 0
      plan%packed_p2p_active = .false.
      plan%ready = .false.
      plan%direct = .false.
      plan%sparse = .false.
      plan%streamed = .false.
      plan%stream_width = default_stream_elements
   end subroutine ! Release_Transfer

   !------------------------------------------------------------------------------
   !! Field movement through the selected transport
   !------------------------------------------------------------------------------
   !! Input and output storage must be disjoint, including same-shape transfers.
   subroutine Move_Real(plan, input, output)
      type(Transfer_Plan), intent(inout) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout), target :: output(:, :, :)
      integer :: rank, ierr, np, n(3), q, i, j, k, xyz(3), idx(3), base, input_shape(3)
      integer(ISFFT_INT64_KIND) :: t, send_size, stride(3), step(3), begin, flat_base
      real(kind=ISFFT_REAL_KIND), pointer :: flat_in(:), flat_out(:)

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      if (.not. plan%ready) error stop 'unbuilt transfer'
      if (any(shape(input) /= Layout_Shape(plan%source, rank))) error stop 'transfer input shape'
      if (plan%chunked .or. plan%streamed) call Bound_Layout_Workspace(min(plan%chunk_limit, plan%stream_width))
      n = Layout_Shape(plan%destination, rank)
      if (allocated(output)) then
         if (any(shape(output) /= n)) deallocate (output)
      end if
      if (.not. allocated(output)) allocate (output(n(1), n(2), n(3)))
      if (Same_Layout(plan%source, plan%destination)) then
         output = input
         return
      end if
      if (plan%local_only) then
         n = shape(output)
         input_shape = shape(input)
         if (size(output, kind=ISFFT_INT64_KIND) == 0) return
         flat_in(1:size(input, kind=ISFFT_INT64_KIND)) => input
         flat_out(1:size(output, kind=ISFFT_INT64_KIND)) => output
         stride(plan%source%axis_order) = &
            [1_ISFFT_INT64_KIND, int(input_shape(1), ISFFT_INT64_KIND), int(input_shape(1), ISFFT_INT64_KIND) * &
               input_shape(2)]
         step = stride(plan%destination%axis_order)
         flat_base = sum(int(plan%send_box_offsets(:, rank), ISFFT_INT64_KIND) * stride)
         t = 0
         do k = 0, n(3) - 1
            do j = 0, n(2) - 1
               begin = flat_base + k * step(3) + j * step(2) + 1
               flat_out(t + 1:t + n(1)) = flat_in(begin:begin + (n(1) - 1) * step(1):step(1))
               t = t + n(1)
            end do
         end do
         return
      end if
      np = size(plan%send_counts)
      base = ISFFT_REAL_MPI_TYPE
      if (plan%packed_p2p_active .and. .not. (plan%chunked .or. plan%streamed) .and. &
         plan%gather_comm == MPI_COMM_NULL) then
         call Packed_P2P_Real(plan, input, output, base)
         return
      end if
      if (plan%sparse .and. .not. (plan%chunked .or. plan%streamed) .and. plan%gather_comm == MPI_COMM_NULL) then
         call Sparse_Real(plan, input, output)
         return
      end if
      if (plan%direct .and. .not. (plan%chunked .or. plan%streamed) .and. &
         .not. plan%local_only .and. plan%gather_comm == MPI_COMM_NULL) then
         call MPI_Alltoallw( &
            input, plan%send_type_counts(plan%peers), plan%send_byte_displs(plan%peers), plan%send_types(plan%peers), &
            output, plan%recv_type_counts(plan%peers), plan%recv_byte_displs(plan%peers), &
            plan%recv_types(plan%peers), plan%exchange_comm, &
            ierr)
         return
      end if
      if (plan%chunked .or. plan%streamed) then
         call Stream_Real(plan, input, output, base)
         return
      end if
      send_size = sum(plan%send_counts)
      if (plan%gather_comm /= MPI_COMM_NULL) send_size = plan%gather_count
      call Reserve_R(send_real_buffer, max(1_ISFFT_INT64_KIND, send_size))
      call Reserve_R(recv_real_buffer, max(1_ISFFT_INT64_KIND, sum(plan%recv_counts)))
      if (plan%gather_comm /= MPI_COMM_NULL) then
         n(plan%source%axis_order) = shape(input)
         t = 0
         !! Serialize one local piece in physical x/y/z order.
         xyz = 0
         do k = 1, n(3)
            do j = 1, n(2)
               do i = 1, n(1)
                  xyz = [i, j, k]
                  idx = xyz(plan%source%axis_order)
                  t = t + 1
                  send_real_buffer(t) = input(idx(1), idx(2), idx(3))
               end do
            end do
         end do
         call Gather_Real(plan, send_real_buffer, recv_real_buffer, base)
      else
         do q = 0, np - 1
            t = plan%send_displs(q)
            do k = 0, plan%send_box_sizes(3, q) - 1
               do j = 0, plan%send_box_sizes(2, q) - 1
                  do i = 0, plan%send_box_sizes(1, q) - 1
                     xyz = plan%send_box_offsets(:, q) + [i, j, k] + 1
                     idx = xyz(plan%source%axis_order)
                     t = t + 1
                     send_real_buffer(t) = input(idx(1), idx(2), idx(3))
                  end do
               end do
            end do
         end do
         if (plan%local_only) then
            recv_real_buffer(1:sum(plan%recv_counts)) = send_real_buffer(1:sum(plan%send_counts))
         else
            call Exchange_Real(plan, send_real_buffer, recv_real_buffer, base)
         end if
      end if
      do q = 0, np - 1
         t = plan%recv_displs(q)
         do k = 0, plan%recv_box_sizes(3, q) - 1
            do j = 0, plan%recv_box_sizes(2, q) - 1
               do i = 0, plan%recv_box_sizes(1, q) - 1
                  xyz = plan%recv_box_offsets(:, q) + [i, j, k] + 1
                  idx = xyz(plan%destination%axis_order)
                  t = t + 1
                  output(idx(1), idx(2), idx(3)) = recv_real_buffer(t)
               end do
            end do
         end do
      end do
   end subroutine ! Move_Real

   !! Keep real and complex transport dispatch, counts and packing order aligned.
   subroutine Move_Complex(plan, input, output)
      type(Transfer_Plan), intent(inout) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable, intent(inout), target :: output(:, :, :)
      integer :: rank, ierr, np, n(3), q, i, j, k, xyz(3), idx(3), base, input_shape(3)
      integer(ISFFT_INT64_KIND) :: t, send_size, stride(3), step(3), begin, flat_base
      complex(kind=ISFFT_REAL_KIND), pointer :: flat_in(:), flat_out(:)

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      if (.not. plan%ready) error stop 'unbuilt transfer'
      if (any(shape(input) /= Layout_Shape(plan%source, rank))) error stop 'transfer input shape'
      if (plan%chunked .or. plan%streamed) call Bound_Layout_Workspace(min(plan%chunk_limit, plan%stream_width))
      n = Layout_Shape(plan%destination, rank)
      if (allocated(output)) then
         if (any(shape(output) /= n)) deallocate (output)
      end if
      if (.not. allocated(output)) allocate (output(n(1), n(2), n(3)))
      if (Same_Layout(plan%source, plan%destination)) then
         output = input
         return
      end if
      if (plan%local_only) then
         n = shape(output)
         input_shape = shape(input)
         if (size(output, kind=ISFFT_INT64_KIND) == 0) return
         flat_in(1:size(input, kind=ISFFT_INT64_KIND)) => input
         flat_out(1:size(output, kind=ISFFT_INT64_KIND)) => output
         stride(plan%source%axis_order) = &
            [1_ISFFT_INT64_KIND, int(input_shape(1), ISFFT_INT64_KIND), int(input_shape(1), ISFFT_INT64_KIND) * &
               input_shape(2)]
         step = stride(plan%destination%axis_order)
         flat_base = sum(int(plan%send_box_offsets(:, rank), ISFFT_INT64_KIND) * stride)
         t = 0
         do k = 0, n(3) - 1
            do j = 0, n(2) - 1
               begin = flat_base + k * step(3) + j * step(2) + 1
               flat_out(t + 1:t + n(1)) = flat_in(begin:begin + (n(1) - 1) * step(1):step(1))
               t = t + n(1)
            end do
         end do
         return
      end if
      np = size(plan%send_counts)
      base = ISFFT_REAL_MPI_TYPE
      call MPI_Type_match_size(MPI_TYPECLASS_COMPLEX, 2 * storage_size(0.0_ISFFT_REAL_KIND) / 8, base, ierr)
      if (plan%packed_p2p_active .and. .not. (plan%chunked .or. plan%streamed) .and. &
         plan%gather_comm == MPI_COMM_NULL) then
         call Packed_P2P_Complex(plan, input, output, base)
         return
      end if
      if (plan%sparse .and. .not. (plan%chunked .or. plan%streamed) .and. plan%gather_comm == MPI_COMM_NULL) then
         call Sparse_Complex(plan, input, output)
         return
      end if
      if (plan%direct .and. .not. (plan%chunked .or. plan%streamed) .and. &
         .not. plan%local_only .and. plan%gather_comm == MPI_COMM_NULL) then
         call MPI_Alltoallw( &
            input, plan%send_type_counts(plan%peers), plan%send_byte_displs(plan%peers), plan%send_types(plan%peers), &
            output, plan%recv_type_counts(plan%peers), plan%recv_byte_displs(plan%peers), &
            plan%recv_types(plan%peers), plan%exchange_comm, &
            ierr)
         return
      end if
      if (plan%chunked .or. plan%streamed) then
         call Stream_Complex(plan, input, output, base)
         return
      end if
      send_size = sum(plan%send_counts)
      if (plan%gather_comm /= MPI_COMM_NULL) send_size = plan%gather_count
      call Reserve_Z(send_complex_buffer, max(1_ISFFT_INT64_KIND, send_size))
      call Reserve_Z(recv_complex_buffer, max(1_ISFFT_INT64_KIND, sum(plan%recv_counts)))
      if (plan%gather_comm /= MPI_COMM_NULL) then
         n(plan%source%axis_order) = shape(input)
         t = 0
         !! Serialize one local piece in physical x/y/z order.
         xyz = 0
         do k = 1, n(3)
            do j = 1, n(2)
               do i = 1, n(1)
                  xyz = [i, j, k]
                  idx = xyz(plan%source%axis_order)
                  t = t + 1
                  send_complex_buffer(t) = input(idx(1), idx(2), idx(3))
               end do
            end do
         end do
         call Gather_Complex(plan, send_complex_buffer, recv_complex_buffer, base)
      else
         do q = 0, np - 1
            t = plan%send_displs(q)
            do k = 0, plan%send_box_sizes(3, q) - 1
               do j = 0, plan%send_box_sizes(2, q) - 1
                  do i = 0, plan%send_box_sizes(1, q) - 1
                     xyz = plan%send_box_offsets(:, q) + [i, j, k] + 1
                     idx = xyz(plan%source%axis_order)
                     t = t + 1
                     send_complex_buffer(t) = input(idx(1), idx(2), idx(3))
                  end do
               end do
            end do
         end do
         if (plan%local_only) then
            recv_complex_buffer(1:sum(plan%recv_counts)) = send_complex_buffer(1:sum(plan%send_counts))
         else
            call Exchange_Complex(plan, send_complex_buffer, recv_complex_buffer, base)
         end if
      end if
      do q = 0, np - 1
         t = plan%recv_displs(q)
         do k = 0, plan%recv_box_sizes(3, q) - 1
            do j = 0, plan%recv_box_sizes(2, q) - 1
               do i = 0, plan%recv_box_sizes(1, q) - 1
                  xyz = plan%recv_box_offsets(:, q) + [i, j, k] + 1
                  idx = xyz(plan%destination%axis_order)
                  t = t + 1
                  output(idx(1), idx(2), idx(3)) = recv_complex_buffer(t)
               end do
            end do
         end do
      end do
   end subroutine ! Move_Complex

   !! Serialize physical x/y/z order in bounded tiles, independently of memory order.
   subroutine Stream_Real(plan, input, output, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      real(kind=ISFFT_REAL_KIND), intent(out), contiguous, target :: output(:, :, :)
      integer, intent(in) :: mpi_datatype
      integer :: rank, np, step, dest, source, ds, ss, ns, nr, ierr, status(MPI_STATUS_SIZE)
      integer :: root, group_rank, c(3), q, nphysical(3), off(3)
      integer(ISFFT_INT64_KIND) :: capacity, offset

      capacity = max(1_ISFFT_INT64_KIND, min(plan%chunk_limit, plan%max_peer, plan%stream_width))
      call Reserve_R(send_real_buffer, capacity)
      call Reserve_R(recv_real_buffer, capacity)
      if (plan%gather_comm /= MPI_COMM_NULL) then
         call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
         call MPI_Comm_rank(plan%gather_comm, group_rank, ierr)
         call MPI_Comm_size(plan%gather_comm, np, ierr)
         call Coordinates(plan%source, rank, c)
         call Layout_Box(plan%source, rank, nphysical, off)
         do root = 0, np - 1
            c(plan%destination%replicated_axis) = root
            q = c(1) + plan%source%process_grid(1) * (c(2) + plan%source%process_grid(2) * c(3))
            do offset = 0_ISFFT_INT64_KIND, plan%gather_counts(root) - 1, capacity
               nr = int(min(capacity, plan%gather_counts(root) - offset))
               if (group_rank == root) &
                  call Pack_Real_Chunk(input, plan%source%axis_order, [0, 0, 0], nphysical, offset, nr, &
                  send_real_buffer)
               call MPI_Bcast(send_real_buffer, nr, mpi_datatype, root, plan%gather_comm, ierr)
               call Unpack_Real_Chunk(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
                  plan%recv_box_sizes(:, q), &
                  offset, nr, send_real_buffer)
            end do
         end do
      else
         call MPI_Comm_rank(plan%exchange_comm, rank, ierr)
         call MPI_Comm_size(plan%exchange_comm, np, ierr)
         do step = 0, np - 1
            dest = mod(rank + step, np)
            source = mod(rank - step + np, np)
            ds = plan%peers(dest + 1)
            ss = plan%peers(source + 1)
            !! Every participant uses the same number of rounds, including zero-count peers.
            do offset = 0_ISFFT_INT64_KIND, plan%max_peer - 1, capacity
               ns = int(min(capacity, max(0_ISFFT_INT64_KIND, plan%send_counts(ds) - offset)))
               nr = int(min(capacity, max(0_ISFFT_INT64_KIND, plan%recv_counts(ss) - offset)))
               if (ns > 0) call Pack_Real_Chunk(input, plan%source%axis_order, plan%send_box_offsets(:, ds), &
                  plan%send_box_sizes(:, ds), offset, ns, send_real_buffer)
               call MPI_Sendrecv(send_real_buffer, ns, mpi_datatype, dest, 0, recv_real_buffer, nr, mpi_datatype, &
                  source, 0, &
                  plan%exchange_comm, status, ierr)
               if (nr > 0) call Unpack_Real_Chunk(output, plan%destination%axis_order, plan%recv_box_offsets(:, ss), &
                  plan%recv_box_sizes(:, ss), offset, nr, recv_real_buffer)
            end do
         end do
      end if
   end subroutine ! Stream_Real

   subroutine Pack_Real_Chunk(field, axis_order, box_offset, box_shape, chunk_offset, element_count, buffer)
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: field(:, :, :)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3), element_count
      integer(ISFFT_INT64_KIND), intent(in) :: chunk_offset
      real(kind=ISFFT_REAL_KIND), intent(out) :: buffer(:)
      real(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: t, width
      integer(ISFFT_INT64_KIND) :: cursor, xyz(3), stride(3), first, last

      if (element_count == 0) return
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      stride(axis_order) = &
         [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
            size(field, 2)]
      t = 1
      cursor = chunk_offset
      do while (t <= element_count)
         xyz(1) = mod(cursor, int(box_shape(1), ISFFT_INT64_KIND))
         xyz(2) = mod(cursor / int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(2), ISFFT_INT64_KIND))
         xyz(3) = cursor / (int(box_shape(1), ISFFT_INT64_KIND) * box_shape(2))
         width = int(min(int(element_count - t + 1, ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) - xyz(1)))
         first = 1 + sum((int(box_offset, ISFFT_INT64_KIND) + xyz) * stride)
         last = first + int(width - 1, ISFFT_INT64_KIND) * stride(1)
         buffer(t:t + width - 1) = flat(first:last:stride(1))
         t = t + width
         cursor = cursor + width
      end do
   end subroutine ! Pack_Real_Chunk

   subroutine Unpack_Real_Chunk(field, axis_order, box_offset, box_shape, chunk_offset, element_count, buffer)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: field(:, :, :)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3), element_count
      integer(ISFFT_INT64_KIND), intent(in) :: chunk_offset
      real(kind=ISFFT_REAL_KIND), intent(in) :: buffer(:)
      real(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: t, width
      integer(ISFFT_INT64_KIND) :: cursor, xyz(3), stride(3), first, last

      if (element_count == 0) return
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      stride(axis_order) = &
         [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
            size(field, 2)]
      t = 1
      cursor = chunk_offset
      do while (t <= element_count)
         xyz(1) = mod(cursor, int(box_shape(1), ISFFT_INT64_KIND))
         xyz(2) = mod(cursor / int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(2), ISFFT_INT64_KIND))
         xyz(3) = cursor / (int(box_shape(1), ISFFT_INT64_KIND) * box_shape(2))
         width = int(min(int(element_count - t + 1, ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) - xyz(1)))
         first = 1 + sum((int(box_offset, ISFFT_INT64_KIND) + xyz) * stride)
         last = first + int(width - 1, ISFFT_INT64_KIND) * stride(1)
         flat(first:last:stride(1)) = buffer(t:t + width - 1)
         t = t + width
         cursor = cursor + width
      end do
   end subroutine ! Unpack_Real_Chunk

   !! Serialize physical x/y/z order in bounded tiles, independently of memory order.
   subroutine Stream_Complex(plan, input, output, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(out), contiguous, target :: output(:, :, :)
      integer, intent(in) :: mpi_datatype
      integer :: rank, np, step, dest, source, ds, ss, ns, nr, ierr, status(MPI_STATUS_SIZE)
      integer :: root, group_rank, c(3), q, nphysical(3), off(3)
      integer(ISFFT_INT64_KIND) :: capacity, offset

      capacity = max(1_ISFFT_INT64_KIND, min(plan%chunk_limit, plan%max_peer, plan%stream_width))
      call Reserve_Z(send_complex_buffer, capacity)
      call Reserve_Z(recv_complex_buffer, capacity)
      if (plan%gather_comm /= MPI_COMM_NULL) then
         call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
         call MPI_Comm_rank(plan%gather_comm, group_rank, ierr)
         call MPI_Comm_size(plan%gather_comm, np, ierr)
         call Coordinates(plan%source, rank, c)
         call Layout_Box(plan%source, rank, nphysical, off)
         do root = 0, np - 1
            c(plan%destination%replicated_axis) = root
            q = c(1) + plan%source%process_grid(1) * (c(2) + plan%source%process_grid(2) * c(3))
            do offset = 0_ISFFT_INT64_KIND, plan%gather_counts(root) - 1, capacity
               nr = int(min(capacity, plan%gather_counts(root) - offset))
               if (group_rank == root) call Pack_Complex_Chunk( &
                  input, plan%source%axis_order, [0, 0, 0], nphysical, offset, nr, send_complex_buffer)
               call MPI_Bcast(send_complex_buffer, nr, mpi_datatype, root, plan%gather_comm, ierr)
               call Unpack_Complex_Chunk(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
                  plan%recv_box_sizes(:, q), offset, nr, send_complex_buffer)
            end do
         end do
      else
         call MPI_Comm_rank(plan%exchange_comm, rank, ierr)
         call MPI_Comm_size(plan%exchange_comm, np, ierr)
         do step = 0, np - 1
            dest = mod(rank + step, np)
            source = mod(rank - step + np, np)
            ds = plan%peers(dest + 1)
            ss = plan%peers(source + 1)
            !! Every participant uses the same number of rounds, including zero-count peers.
            do offset = 0_ISFFT_INT64_KIND, plan%max_peer - 1, capacity
               ns = int(min(capacity, max(0_ISFFT_INT64_KIND, plan%send_counts(ds) - offset)))
               nr = int(min(capacity, max(0_ISFFT_INT64_KIND, plan%recv_counts(ss) - offset)))
               if (ns > 0) call Pack_Complex_Chunk(input, plan%source%axis_order, plan%send_box_offsets(:, ds), &
                  plan%send_box_sizes(:, ds), offset, ns, send_complex_buffer)
               call MPI_Sendrecv(send_complex_buffer, ns, mpi_datatype, dest, 0, recv_complex_buffer, nr, &
                  mpi_datatype, source, 0, &
                  plan%exchange_comm, status, ierr)
               if (nr > 0) call Unpack_Complex_Chunk( &
                  output, plan%destination%axis_order, plan%recv_box_offsets(:, ss), plan%recv_box_sizes(:, ss), &
                  offset, nr, &
                  recv_complex_buffer)
            end do
         end do
      end if
   end subroutine ! Stream_Complex

   subroutine Pack_Complex_Chunk(field, axis_order, box_offset, box_shape, chunk_offset, element_count, buffer)
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: field(:, :, :)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3), element_count
      integer(ISFFT_INT64_KIND), intent(in) :: chunk_offset
      complex(kind=ISFFT_REAL_KIND), intent(out) :: buffer(:)
      complex(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: t, width
      integer(ISFFT_INT64_KIND) :: cursor, xyz(3), stride(3), first, last

      if (element_count == 0) return
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      stride(axis_order) = &
         [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
            size(field, 2)]
      t = 1
      cursor = chunk_offset
      do while (t <= element_count)
         xyz(1) = mod(cursor, int(box_shape(1), ISFFT_INT64_KIND))
         xyz(2) = mod(cursor / int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(2), ISFFT_INT64_KIND))
         xyz(3) = cursor / (int(box_shape(1), ISFFT_INT64_KIND) * box_shape(2))
         width = int(min(int(element_count - t + 1, ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) - xyz(1)))
         first = 1 + sum((int(box_offset, ISFFT_INT64_KIND) + xyz) * stride)
         last = first + int(width - 1, ISFFT_INT64_KIND) * stride(1)
         buffer(t:t + width - 1) = flat(first:last:stride(1))
         t = t + width
         cursor = cursor + width
      end do
   end subroutine ! Pack_Complex_Chunk

   subroutine Unpack_Complex_Chunk(field, axis_order, box_offset, box_shape, chunk_offset, element_count, buffer)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: field(:, :, :)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3), element_count
      integer(ISFFT_INT64_KIND), intent(in) :: chunk_offset
      complex(kind=ISFFT_REAL_KIND), intent(in) :: buffer(:)
      complex(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: t, width
      integer(ISFFT_INT64_KIND) :: cursor, xyz(3), stride(3), first, last

      if (element_count == 0) return
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      stride(axis_order) = &
         [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
            size(field, 2)]
      t = 1
      cursor = chunk_offset
      do while (t <= element_count)
         xyz(1) = mod(cursor, int(box_shape(1), ISFFT_INT64_KIND))
         xyz(2) = mod(cursor / int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(2), ISFFT_INT64_KIND))
         xyz(3) = cursor / (int(box_shape(1), ISFFT_INT64_KIND) * box_shape(2))
         width = int(min(int(element_count - t + 1, ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) - xyz(1)))
         first = 1 + sum((int(box_offset, ISFFT_INT64_KIND) + xyz) * stride)
         last = first + int(width - 1, ISFFT_INT64_KIND) * stride(1)
         flat(first:last:stride(1)) = buffer(t:t + width - 1)
         t = t + width
         cursor = cursor + width
      end do
   end subroutine ! Unpack_Complex_Chunk

   !------------------------------------------------------------------------------
   !! Stable transport method dispatch
   !------------------------------------------------------------------------------
   !! The chunked flag enforces MPI count safety and takes priority over a method choice.
   !! Stable method IDs: packed_a2a=0, datatype_a2a=1, packed_p2p=2, datatype_p2p=3, memorysave_p2p=4.
   !! Legacy direct/sparse/streamed flags remain compatible with existing diagnostic consumers.
   !! Count-safety and replication constraints take precedence over any requested method.
   subroutine Select_Transfer_Method(plan, method)
      type(Transfer_Plan), intent(inout) :: plan
      integer, intent(in) :: method

      plan%packed_p2p_active = .false.
      plan%direct = .false.
      plan%sparse = .false.
      plan%streamed = plan%chunked
      if (plan%local_only .or. plan%chunked) return
      select case (method)
      case (packed_a2a)
      case (datatype_a2a)
         plan%direct = plan%gather_comm == MPI_COMM_NULL
      case (memorysave_p2p)
         plan%streamed = .true.
      case (datatype_p2p)
         plan%sparse = plan%gather_comm == MPI_COMM_NULL .and. allocated(plan%requests)
      case (packed_p2p)
         plan%packed_p2p_active = plan%gather_comm == MPI_COMM_NULL .and. allocated(plan%requests)
      case default
         error stop 'invalid transfer method'
      end select
   end subroutine ! Select_Transfer_Method

   function Transfer_Method_Name(method) result(name)
      integer, intent(in) :: method
      character(24) :: name

      select case (method)
      case (packed_a2a)
         name = 'packed_a2a'
      case (datatype_a2a)
         name = 'datatype_a2a'
      case (memorysave_p2p)
         name = 'memorysave_p2p'
      case (datatype_p2p)
         name = 'datatype_p2p'
      case (packed_p2p)
         name = 'packed_p2p'
      case default
         name = 'unknown'
      end select
   end function ! Transfer_Method_Name

   !------------------------------------------------------------------------------
   !! Packed nonblocking point-to-point transport
   !------------------------------------------------------------------------------
   !! Receives are posted before packing. Each send owns a disjoint buffer slice;
   !! completed receives are unpacked while other messages may still be in flight.
   !! No execution-time request allocation; all sends complete before scratch reuse.
   subroutine Packed_P2P_Real(plan, input, output, mpi_datatype)
      type(Transfer_Plan), intent(inout) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: output(:, :, :)
      integer, intent(in) :: mpi_datatype
      integer :: i, peer, q, nrecv, nsend, done, ierr, rank, nremote
      integer(ISFFT_INT64_KIND) :: first, last

      call Reserve_R(send_real_buffer, max(1_ISFFT_INT64_KIND, sum(plan%send_counts)))
      call Reserve_R(recv_real_buffer, max(1_ISFFT_INT64_KIND, sum(plan%recv_counts)))
      nrecv = size(plan%sparse_recv)
      nsend = size(plan%sparse_send)
      call MPI_Comm_rank(plan%exchange_comm, rank, ierr)
      plan%requests = MPI_REQUEST_NULL
      nremote = 0

      !! 1. Post remote receives
      do i = 1, nrecv
         peer = plan%sparse_recv(i)
         q = plan%peers(peer + 1)
         if (peer == rank) cycle
         nremote = nremote + 1
         first = plan%recv_displs(q) + 1
         last = plan%recv_displs(q) + plan%recv_counts(q)
         call MPI_Irecv(recv_real_buffer(first:last), int(plan%recv_counts(q)), mpi_datatype, peer, 0, &
            plan%exchange_comm, &
            plan%requests(i), ierr)
      end do

      !! 2. Pack disjoint send slices and handle self locally
      do i = 1, nsend
         peer = plan%sparse_send(i)
         q = plan%peers(peer + 1)
         first = plan%send_displs(q) + 1
         last = plan%send_displs(q) + plan%send_counts(q)
         call Pack_P2P_Real_Box(input, plan%source%axis_order, plan%send_box_offsets(:, q), plan%send_box_sizes(:, q), &
            send_real_buffer(first:last))
         if (peer == rank) then
            call Unpack_P2P_Real_Box(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
               plan%recv_box_sizes(:, q), &
               send_real_buffer(first:last))
            cycle
         end if
         call MPI_Isend(send_real_buffer(first:last), int(plan%send_counts(q)), mpi_datatype, peer, 0, &
            plan%exchange_comm, &
            plan%requests(nrecv + i), ierr)
      end do

      !! 3. Unpack completed receives
      do i = 1, nremote
         call MPI_Waitany(nrecv, plan%requests(:nrecv), done, MPI_STATUS_IGNORE, ierr)
         if (done == MPI_UNDEFINED) error stop 'packed_p2p incomplete receive'
         peer = plan%sparse_recv(done)
         q = plan%peers(peer + 1)
         first = plan%recv_displs(q) + 1
         last = plan%recv_displs(q) + plan%recv_counts(q)
         call Unpack_P2P_Real_Box(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
            plan%recv_box_sizes(:, q), &
            recv_real_buffer(first:last))
      end do

      !! 4. Complete sends before reusing shared scratch
      call MPI_Waitall(nsend, plan%requests(nrecv + 1:nrecv + nsend), MPI_STATUSES_IGNORE, ierr)
   end subroutine ! Packed_P2P_Real

   !! Receives are posted before packing. Each send owns a disjoint buffer slice;
   !! completed receives are unpacked while other messages may still be in flight.
   !! No execution-time request allocation; all sends complete before scratch reuse.
   subroutine Packed_P2P_Complex(plan, input, output, mpi_datatype)
      type(Transfer_Plan), intent(inout) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: input(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: output(:, :, :)
      integer, intent(in) :: mpi_datatype
      integer :: i, peer, q, nrecv, nsend, done, ierr, rank, nremote
      integer(ISFFT_INT64_KIND) :: first, last

      call Reserve_Z(send_complex_buffer, max(1_ISFFT_INT64_KIND, sum(plan%send_counts)))
      call Reserve_Z(recv_complex_buffer, max(1_ISFFT_INT64_KIND, sum(plan%recv_counts)))
      nrecv = size(plan%sparse_recv)
      nsend = size(plan%sparse_send)
      call MPI_Comm_rank(plan%exchange_comm, rank, ierr)
      plan%requests = MPI_REQUEST_NULL
      nremote = 0

      !! 1. Post remote receives
      do i = 1, nrecv
         peer = plan%sparse_recv(i)
         q = plan%peers(peer + 1)
         if (peer == rank) cycle
         nremote = nremote + 1
         first = plan%recv_displs(q) + 1
         last = plan%recv_displs(q) + plan%recv_counts(q)
         call MPI_Irecv(recv_complex_buffer(first:last), int(plan%recv_counts(q)), mpi_datatype, peer, 0, &
            plan%exchange_comm, &
            plan%requests(i), ierr)
      end do

      !! 2. Pack disjoint send slices and handle self locally
      do i = 1, nsend
         peer = plan%sparse_send(i)
         q = plan%peers(peer + 1)
         first = plan%send_displs(q) + 1
         last = plan%send_displs(q) + plan%send_counts(q)
         call Pack_P2P_Complex_Box(input, plan%source%axis_order, plan%send_box_offsets(:, q), &
            plan%send_box_sizes(:, q), &
            send_complex_buffer(first:last))
         if (peer == rank) then
            call Unpack_P2P_Complex_Box(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
               plan%recv_box_sizes(:, q), send_complex_buffer(first:last))
            cycle
         end if
         call MPI_Isend(send_complex_buffer(first:last), int(plan%send_counts(q)), mpi_datatype, peer, 0, &
            plan%exchange_comm, &
            plan%requests(nrecv + i), ierr)
      end do

      !! 3. Unpack completed receives
      do i = 1, nremote
         call MPI_Waitany(nrecv, plan%requests(:nrecv), done, MPI_STATUS_IGNORE, ierr)
         if (done == MPI_UNDEFINED) error stop 'packed_p2p incomplete receive'
         peer = plan%sparse_recv(done)
         q = plan%peers(peer + 1)
         first = plan%recv_displs(q) + 1
         last = plan%recv_displs(q) + plan%recv_counts(q)
         call Unpack_P2P_Complex_Box(output, plan%destination%axis_order, plan%recv_box_offsets(:, q), &
            plan%recv_box_sizes(:, q), &
            recv_complex_buffer(first:last))
      end do

      !! 4. Complete sends before reusing shared scratch
      call MPI_Waitall(nsend, plan%requests(nrecv + 1:nrecv + nsend), MPI_STATUSES_IGNORE, ierr)
   end subroutine ! Packed_P2P_Complex

   !------------------------------------------------------------------------------
   !! Storage-order-aware point-to-point packing
   !------------------------------------------------------------------------------
   !! Tile differing storage orders so a strided gather/scatter does not sweep
   !! an entire large plane between neighboring cache-line accesses.
   subroutine Pack_P2P_Real_Box(field, axis_order, box_offset, box_shape, buffer)
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: field(:, :, :)
      real(kind=ISFFT_REAL_KIND), intent(out) :: buffer(:)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3)
      real(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: dims(3), ib, jb, kb, j, k, width
      integer(ISFFT_INT64_KIND) :: stride(3), step(3), offset, start, first, last

      if (size(buffer) == 0) return
      if (axis_order(1) == 1) then
         call Pack_Real_Chunk(field, axis_order, box_offset, box_shape, 0_ISFFT_INT64_KIND, size(buffer), buffer)
         return
      end if
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      dims = box_shape(axis_order)
      stride = [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
         size(field, 2)]
      step = [1_ISFFT_INT64_KIND, int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) * &
         box_shape(2)]
      step = step(axis_order)
      offset = sum(int(box_offset(axis_order), ISFFT_INT64_KIND) * stride)
      do kb = 0, dims(3) - 1, packing_tile_width
         do jb = 0, dims(2) - 1, packing_tile_width
            do ib = 0, dims(1) - 1, packing_tile_width
               width = min(packing_tile_width, dims(1) - ib)
               do k = kb, min(kb + (packing_tile_width - 1), dims(3) - 1)
                  do j = jb, min(jb + (packing_tile_width - 1), dims(2) - 1)
                     start = 1 + offset + ib + int(j, ISFFT_INT64_KIND) * stride(2) + int(k, ISFFT_INT64_KIND) * &
                        stride(3)
                     first = 1 + int(ib, ISFFT_INT64_KIND) * step(1) + int(j, ISFFT_INT64_KIND) * step(2) + &
                        int(k, ISFFT_INT64_KIND) * step(3)
                     last = first + int(width - 1, ISFFT_INT64_KIND) * step(1)
                     buffer(first:last:step(1)) = flat(start:start + width - 1)
                  end do
               end do
            end do
         end do
      end do
   end subroutine ! Pack_P2P_Real_Box

   !! Tile differing storage orders so a strided gather/scatter does not sweep
   !! an entire large plane between neighboring cache-line accesses.
   subroutine Unpack_P2P_Real_Box(field, axis_order, box_offset, box_shape, buffer)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: field(:, :, :)
      real(kind=ISFFT_REAL_KIND), intent(in) :: buffer(:)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3)
      real(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: dims(3), ib, jb, kb, j, k, width
      integer(ISFFT_INT64_KIND) :: stride(3), step(3), offset, start, first, last

      if (size(buffer) == 0) return
      if (axis_order(1) == 1) then
         call Unpack_Real_Chunk(field, axis_order, box_offset, box_shape, 0_ISFFT_INT64_KIND, size(buffer), buffer)
         return
      end if
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      dims = box_shape(axis_order)
      stride = [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
         size(field, 2)]
      step = [1_ISFFT_INT64_KIND, int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) * &
         box_shape(2)]
      step = step(axis_order)
      offset = sum(int(box_offset(axis_order), ISFFT_INT64_KIND) * stride)
      do kb = 0, dims(3) - 1, packing_tile_width
         do jb = 0, dims(2) - 1, packing_tile_width
            do ib = 0, dims(1) - 1, packing_tile_width
               width = min(packing_tile_width, dims(1) - ib)
               do k = kb, min(kb + (packing_tile_width - 1), dims(3) - 1)
                  do j = jb, min(jb + (packing_tile_width - 1), dims(2) - 1)
                     start = 1 + offset + ib + int(j, ISFFT_INT64_KIND) * stride(2) + int(k, ISFFT_INT64_KIND) * &
                        stride(3)
                     first = 1 + int(ib, ISFFT_INT64_KIND) * step(1) + int(j, ISFFT_INT64_KIND) * step(2) + &
                        int(k, ISFFT_INT64_KIND) * step(3)
                     last = first + int(width - 1, ISFFT_INT64_KIND) * step(1)
                     flat(start:start + width - 1) = buffer(first:last:step(1))
                  end do
               end do
            end do
         end do
      end do
   end subroutine ! Unpack_P2P_Real_Box

   !! Tile differing storage orders so a strided gather/scatter does not sweep
   !! an entire large plane between neighboring cache-line accesses.
   subroutine Pack_P2P_Complex_Box(field, axis_order, box_offset, box_shape, buffer)
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: field(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(out) :: buffer(:)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3)
      complex(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: dims(3), ib, jb, kb, j, k, width
      integer(ISFFT_INT64_KIND) :: stride(3), step(3), offset, start, first, last

      if (size(buffer) == 0) return
      if (axis_order(1) == 1) then
         call Pack_Complex_Chunk(field, axis_order, box_offset, box_shape, 0_ISFFT_INT64_KIND, size(buffer), buffer)
         return
      end if
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      dims = box_shape(axis_order)
      stride = [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
         size(field, 2)]
      step = [1_ISFFT_INT64_KIND, int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) * &
         box_shape(2)]
      step = step(axis_order)
      offset = sum(int(box_offset(axis_order), ISFFT_INT64_KIND) * stride)
      do kb = 0, dims(3) - 1, packing_tile_width
         do jb = 0, dims(2) - 1, packing_tile_width
            do ib = 0, dims(1) - 1, packing_tile_width
               width = min(packing_tile_width, dims(1) - ib)
               do k = kb, min(kb + (packing_tile_width - 1), dims(3) - 1)
                  do j = jb, min(jb + (packing_tile_width - 1), dims(2) - 1)
                     start = 1 + offset + ib + int(j, ISFFT_INT64_KIND) * stride(2) + int(k, ISFFT_INT64_KIND) * &
                        stride(3)
                     first = 1 + int(ib, ISFFT_INT64_KIND) * step(1) + int(j, ISFFT_INT64_KIND) * step(2) + &
                        int(k, ISFFT_INT64_KIND) * step(3)
                     last = first + int(width - 1, ISFFT_INT64_KIND) * step(1)
                     buffer(first:last:step(1)) = flat(start:start + width - 1)
                  end do
               end do
            end do
         end do
      end do
   end subroutine ! Pack_P2P_Complex_Box

   !! Tile differing storage orders so a strided gather/scatter does not sweep
   !! an entire large plane between neighboring cache-line accesses.
   subroutine Unpack_P2P_Complex_Box(field, axis_order, box_offset, box_shape, buffer)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: field(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(in) :: buffer(:)
      integer, intent(in) :: axis_order(3), box_offset(3), box_shape(3)
      complex(kind=ISFFT_REAL_KIND), pointer :: flat(:)
      integer :: dims(3), ib, jb, kb, j, k, width
      integer(ISFFT_INT64_KIND) :: stride(3), step(3), offset, start, first, last

      if (size(buffer) == 0) return
      if (axis_order(1) == 1) then
         call Unpack_Complex_Chunk(field, axis_order, box_offset, box_shape, 0_ISFFT_INT64_KIND, size(buffer), buffer)
         return
      end if
      flat(1:size(field, kind=ISFFT_INT64_KIND)) => field
      dims = box_shape(axis_order)
      stride = [1_ISFFT_INT64_KIND, int(size(field, 1), ISFFT_INT64_KIND), int(size(field, 1), ISFFT_INT64_KIND) * &
         size(field, 2)]
      step = [1_ISFFT_INT64_KIND, int(box_shape(1), ISFFT_INT64_KIND), int(box_shape(1), ISFFT_INT64_KIND) * &
         box_shape(2)]
      step = step(axis_order)
      offset = sum(int(box_offset(axis_order), ISFFT_INT64_KIND) * stride)
      do kb = 0, dims(3) - 1, packing_tile_width
         do jb = 0, dims(2) - 1, packing_tile_width
            do ib = 0, dims(1) - 1, packing_tile_width
               width = min(packing_tile_width, dims(1) - ib)
               do k = kb, min(kb + (packing_tile_width - 1), dims(3) - 1)
                  do j = jb, min(jb + (packing_tile_width - 1), dims(2) - 1)
                     start = 1 + offset + ib + int(j, ISFFT_INT64_KIND) * stride(2) + int(k, ISFFT_INT64_KIND) * &
                        stride(3)
                     first = 1 + int(ib, ISFFT_INT64_KIND) * step(1) + int(j, ISFFT_INT64_KIND) * step(2) + &
                        int(k, ISFFT_INT64_KIND) * step(3)
                     last = first + int(width - 1, ISFFT_INT64_KIND) * step(1)
                     flat(start:start + width - 1) = buffer(first:last:step(1))
                  end do
               end do
            end do
         end do
      end do
   end subroutine ! Unpack_P2P_Complex_Box

   !------------------------------------------------------------------------------
   !! Derived-datatype point-to-point transport
   !------------------------------------------------------------------------------
   subroutine Sparse_Real(plan, input, output)
      type(Transfer_Plan), intent(inout) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, asynchronous :: input(:, :, :)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, asynchronous :: output(:, :, :)
      integer :: i, peer, q, nreq, ierr

      nreq = 0
      do i = 1, size(plan%sparse_recv)
         peer = plan%sparse_recv(i)
         q = plan%peers(peer + 1)
         nreq = nreq + 1
         call MPI_Irecv(output, 1, plan%recv_types(q), peer, 0, plan%exchange_comm, plan%requests(nreq), ierr)
      end do
      do i = 1, size(plan%sparse_send)
         peer = plan%sparse_send(i)
         q = plan%peers(peer + 1)
         nreq = nreq + 1
         call MPI_Isend(input, 1, plan%send_types(q), peer, 0, plan%exchange_comm, plan%requests(nreq), ierr)
      end do
      call MPI_Waitall(nreq, plan%requests, MPI_STATUSES_IGNORE, ierr)
   end subroutine ! Sparse_Real

   subroutine Sparse_Complex(plan, input, output)
      type(Transfer_Plan), intent(inout) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous, asynchronous :: input(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, asynchronous :: output(:, :, :)
      integer :: i, peer, q, nreq, ierr

      nreq = 0
      do i = 1, size(plan%sparse_recv)
         peer = plan%sparse_recv(i)
         q = plan%peers(peer + 1)
         nreq = nreq + 1
         call MPI_Irecv(output, 1, plan%recv_types(q), peer, 0, plan%exchange_comm, plan%requests(nreq), ierr)
      end do
      do i = 1, size(plan%sparse_send)
         peer = plan%sparse_send(i)
         q = plan%peers(peer + 1)
         nreq = nreq + 1
         call MPI_Isend(input, 1, plan%send_types(q), peer, 0, plan%exchange_comm, plan%requests(nreq), ierr)
      end do
      call MPI_Waitall(nreq, plan%requests, MPI_STATUSES_IGNORE, ierr)
   end subroutine ! Sparse_Complex

   !------------------------------------------------------------------------------
   !! Packed collective exchange helpers
   !------------------------------------------------------------------------------
   !! Ordinary-size transfers retain native MPI collectives.
   subroutine Exchange_Real(plan, send, recv, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous :: send(:)
      real(kind=ISFFT_REAL_KIND), intent(out), contiguous :: recv(:)
      integer, intent(in) :: mpi_datatype
      integer :: ierr

      call MPI_Alltoallv( &
         send, int(plan%send_counts(plan%peers)), int(plan%send_displs(plan%peers)), mpi_datatype, &
         recv, int(plan%recv_counts(plan%peers)), int(plan%recv_displs(plan%peers)), mpi_datatype, &
         plan%exchange_comm, ierr)
   end subroutine ! Exchange_Real

   subroutine Gather_Real(plan, send, recv, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous :: send(:)
      real(kind=ISFFT_REAL_KIND), intent(out), contiguous :: recv(:)
      integer, intent(in) :: mpi_datatype
      integer :: ierr

      call MPI_Allgatherv(send, int(plan%gather_count), mpi_datatype, recv, int(plan%gather_counts), &
         int(plan%gather_displs), mpi_datatype, &
         plan%gather_comm, ierr)
   end subroutine ! Gather_Real

   subroutine Exchange_Complex(plan, send, recv, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous :: send(:)
      complex(kind=ISFFT_REAL_KIND), intent(out), contiguous :: recv(:)
      integer, intent(in) :: mpi_datatype
      integer :: ierr

      call MPI_Alltoallv( &
         send, int(plan%send_counts(plan%peers)), int(plan%send_displs(plan%peers)), mpi_datatype, &
         recv, int(plan%recv_counts(plan%peers)), int(plan%recv_displs(plan%peers)), mpi_datatype, &
         plan%exchange_comm, ierr)
   end subroutine ! Exchange_Complex

   subroutine Gather_Complex(plan, send, recv, mpi_datatype)
      type(Transfer_Plan), intent(in) :: plan
      complex(kind=ISFFT_REAL_KIND), intent(in), contiguous :: send(:)
      complex(kind=ISFFT_REAL_KIND), intent(out), contiguous :: recv(:)
      integer, intent(in) :: mpi_datatype
      integer :: ierr

      call MPI_Allgatherv(send, int(plan%gather_count), mpi_datatype, recv, int(plan%gather_counts), &
         int(plan%gather_displs), mpi_datatype, &
         plan%gather_comm, ierr)
   end subroutine ! Gather_Complex

   !! Complete the transfer before parking the source and adopting destination storage.
   subroutine Move_State(plan, field)
      type(Transfer_Plan), intent(inout) :: plan
      type(Spectral_Field), intent(inout) :: field
      real(kind=ISFFT_REAL_KIND), allocatable :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: z(:, :, :)
      integer :: rank, ierr, n(3)

      if (Same_Layout(plan%source, plan%destination)) return
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      n = Layout_Shape(plan%destination, rank)
      if (plan%source%complex_data) then
         call Get_Complex(field, n, z)
         call Move_Field(plan, field%z, z)
         call Park_Field(field)
         call move_alloc(z, field%z)
      else
         call Get_Real(field, n, r)
         call Move_Field(plan, field%r, r)
         call Park_Field(field)
         call move_alloc(r, field%r)
      end if
   end subroutine ! Move_State

   !! Check communication plans before any identity-layout early return.
   subroutine Require_Redistribution_Initialization()

      if (.not. redistribution_forward_plan%ready .or. .not. redistribution_inverse_plan%ready) &
         error stop 'ISFFT redistribution: call ISFFT_Initialize before redistribution or layout activation'
   end subroutine ! Require_Redistribution_Initialization

   !------------------------------------------------------------------------------
   !! Field movement and allocation ownership
   !------------------------------------------------------------------------------
   !! Move data using the old topology before activating destination sizes and offsets.
   !! move_alloc transfers ownership. Scratch and caller storage must stay disjoint.
   subroutine Redistribute_Field(field, to_internal)
      use ISFFT_MPI_Part, only: ISFFT_part_change
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: field(:, :, :)
      logical, intent(in) :: to_internal
      real(kind=ISFFT_REAL_KIND), allocatable :: moved(:, :, :), with_halo(:, :, :), owned(:, :, :)
      integer :: grid(3), halo, n(3)

      call Require_Redistribution_Initialization()
      grid = original_layout%process_grid
      if (to_internal) grid = internal_layout%process_grid
      if (all([npx0, npy0, npz0] == grid)) return
      halo = 0
      if (halo_width > 0 .and. all(lbound(field) == 1 - halo_width) .and. &
         all(shape(field) == [nx, ny, nz] + 2 * halo_width)) halo = halo_width
      if (any(shape(field) /= [nx, ny, nz] + 2 * halo)) error stop 'redistribution input shape'
      n = Layout_Shape(redistribution_inverse_plan%destination, my_id)
      if (to_internal) n = Layout_Shape(redistribution_forward_plan%destination, my_id)

      !! 1. Move into destination storage using the current topology
      call Get_Real(redistribution_field, n, moved)
      if (halo > 0) then
         call Get_Real(redistribution_field, [nx, ny, nz], owned)
         owned = field(1:nx, 1:ny, 1:nz)
         if (to_internal) then
            call Move_Field(redistribution_forward_plan, owned, moved)
         else
            call Move_Field(redistribution_inverse_plan, owned, moved)
         end if
         call move_alloc(owned, redistribution_field%r)
         call Park_Field(redistribution_field)
      else
         if (to_internal) then
            call Move_Field(redistribution_forward_plan, field, moved)
         else
            call Move_Field(redistribution_inverse_plan, field, moved)
         end if
      end if

      !! 2. Activate destination topology after field movement
      call ISFFT_part_change(grid(1), grid(2), grid(3))

      !! 3. Transfer allocation ownership and restore halo storage
      if (halo > 0) then
         call Get_Real(redistribution_field, [nx, ny, nz] + 2 * halo, with_halo, 1 - halo)
         with_halo = 0.0_ISFFT_REAL_KIND
         with_halo(1:nx, 1:ny, 1:nz) = moved
         call move_alloc(field, redistribution_field%r)
         call Park_Field(redistribution_field)
         call move_alloc(with_halo, field)
         call move_alloc(moved, redistribution_field%r)
         call Park_Field(redistribution_field)
      else
         call move_alloc(field, redistribution_field%r)
         call Park_Field(redistribution_field)
         call move_alloc(moved, field)
      end if
   end subroutine Redistribute_Field

   subroutine Release_Transform_Plan(plan)
      type(Transform_Plan), intent(inout) :: plan
      integer :: a

      do a = 1, 3
         plan%axes(a)%tail_kernel = Kernel_Plan()
         plan%axes(a)%tiled_block = .false.
         plan%axes(a)%communication_batch_width = 1
         call Release_Transfer(plan%axes(a)%pack)
         call Release_Transfer(plan%axes(a)%unpack)
         call Release_Transfer(plan%axes(a)%reverse_pack)
         call Release_Transfer(plan%axes(a)%reverse_unpack)
         call Release_Transfer(plan%axes(a)%before)
         call Release_Transfer(plan%axes(a)%after)
         call Release_Transfer(plan%axes(a)%reverse_before)
         call Release_Transfer(plan%axes(a)%reverse_after)
      end do
      plan = Transform_Plan()
   end subroutine ! Release_Transform_Plan

   subroutine Begin_Block_Transfer(plan, field, sign, transfer)
      type(Axis_Plan), intent(in) :: plan
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: sign
      type(Block_Transfer), intent(out) :: transfer
      integer :: rank, ierr, a, q

      associate (comm => transfer%comm, parts => transfer%parts, axis => transfer%axis, &
         source => transfer%source, destination => transfer%destination, nin => transfer%nin, &
         nout => transfer%nout, oi => transfer%oi, oo => transfer%oo, &
         others => transfer%others, lines => transfer%lines, capacity => transfer%capacity, &
         local_n => transfer%local_n, full_n => transfer%full_n, out_n => transfer%out_n, &
         result => transfer%result)
         comm = plan%before%gather_comm
         axis = plan%physical_axis
         source = plan%before%source
         destination = plan%after%destination
         if (sign < 0) then
            source = plan%after%destination
            destination = plan%before%source
         end if
         call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
         call MPI_Comm_size(comm, parts, ierr)
         call Layout_Box(source, rank, nin, oi)
         call Layout_Box(destination, rank, nout, oo)
         a = 0
         do q = 1, 3
            if (q == axis) cycle
            a = a + 1
            others(a) = q
         end do
         lines = nin(others(1)) * nin(others(2))
         capacity = max(1, min(plan%communication_batch_width, max(1, lines)))
         local_n = nin(axis)
         full_n = source%global_shape(axis)
         out_n = destination%global_shape(axis)
         if (destination%complex_data) then
            call Get_Complex(field, nout, result%z)
         else
            call Get_Real(field, nout, result%r)
         end if
         if (source%complex_data) then
            transfer%sz => block_send_z
            transfer%rz => block_recv_z
            transfer%tile_z(1:full_n, 1:capacity, 1:1) => block_tile_z(1:full_n * capacity)
         else
            transfer%sr => block_send_r
            transfer%rr => block_recv_r
            transfer%tile_r(1:full_n, 1:capacity, 1:1) => block_tile_r(1:full_n * capacity)
         end if
         if (plan%transform_kind == 1) then
            if (sign > 0) then
               transfer%tile_z(1:out_n, 1:capacity, 1:1) => block_tile_z(1:out_n * capacity)
            else
               transfer%tile_r(1:out_n, 1:capacity, 1:1) => block_tile_r(1:out_n * capacity)
            end if
         end if
      end associate
   end subroutine Begin_Block_Transfer

   subroutine Gather_Block_Tile(transfer, field, first, width)
      type(Block_Transfer), intent(inout) :: transfer
      type(Spectral_Field), intent(in) :: field
      integer, intent(in) :: first, width
      integer :: q, count, line, id, xyz(3), i, position, ierr, partoff, start

      associate (comm => transfer%comm, parts => transfer%parts, axis => transfer%axis, &
         source => transfer%source, nin => transfer%nin, others => transfer%others, &
         local_n => transfer%local_n, full_n => transfer%full_n, sr => transfer%sr, &
         rr => transfer%rr, sz => transfer%sz, rz => transfer%rz, &
         tile_r => transfer%tile_r, tile_z => transfer%tile_z)
            do q = 0, parts - 1
               count = full_n / parts
               if (q < mod(full_n, parts)) count = count + 1
               block_counts(q) = count * width
               block_displs(q) = (q * (full_n / parts) + min(q, mod(full_n, parts))) * width
            end do

            !! Pack local tile lines
            do line = 1, width
               id = first + line - 1
               xyz = 1
               xyz(others(1)) = mod(id, nin(others(1))) + 1
               xyz(others(2)) = id / nin(others(1)) + 1
               do i = 1, local_n
                  xyz(axis) = i
                  position = (line - 1) * local_n + i
                  if (source%complex_data) then
                     sz(position) = field%z(xyz(1), xyz(2), xyz(3))
                  else
                     sr(position) = field%r(xyz(1), xyz(2), xyz(3))
                  end if
               end do
            end do

            !! Gather each rank contribution into complete transform lines
            if (source%complex_data) then
               call MPI_Allgatherv(sz, local_n * width, ISFFT_COMPLEX_MPI_TYPE, rz, block_counts, block_displs, &
                  ISFFT_COMPLEX_MPI_TYPE, comm, ierr)
            else
               call MPI_Allgatherv(sr, local_n * width, ISFFT_REAL_MPI_TYPE, rr, block_counts, block_displs, &
                  ISFFT_REAL_MPI_TYPE, comm, ierr)
            end if
            do q = 0, parts - 1
               count = block_counts(q) / width
               partoff = q * (full_n / parts) + min(q, mod(full_n, parts))
               do line = 1, width
                  start = block_displs(q) + (line - 1) * count
                  if (source%complex_data) then
                     tile_z(partoff + 1:partoff + count, line, 1) = rz(start + 1:start + count)
                  else
                     tile_r(partoff + 1:partoff + count, line, 1) = rr(start + 1:start + count)
                  end if
               end do
            end do

      end associate
   end subroutine Gather_Block_Tile

   subroutine Store_Block_Tile(transfer, first, width)
      type(Block_Transfer), intent(inout) :: transfer
      integer, intent(in) :: first, width
      integer :: line, id, xyz(3), i

      associate (axis => transfer%axis, source => transfer%source, destination => transfer%destination, &
         nin => transfer%nin, nout => transfer%nout, oo => transfer%oo, &
         others => transfer%others, result => transfer%result, tile_r => transfer%tile_r, &
         tile_z => transfer%tile_z)
            !! Scatter only the destination-owned portion
            do line = 1, width
               id = first + line - 1
               xyz = 1
               xyz(others(1)) = mod(id, nin(others(1))) + 1
               xyz(others(2)) = id / nin(others(1)) + 1
               do i = 1, nout(axis)
                  xyz(axis) = i
                  if (destination%complex_data) then
                     result%z(xyz(1), xyz(2), xyz(3)) = tile_z(oo(axis) + i, line, 1)
                  else
                     result%r(xyz(1), xyz(2), xyz(3)) = tile_r(oo(axis) + i, line, 1)
                  end if
               end do
            end do
      end associate
   end subroutine Store_Block_Tile

   subroutine End_Block_Transfer(transfer, field)
      type(Block_Transfer), intent(inout) :: transfer
      type(Spectral_Field), intent(inout) :: field

      associate (result => transfer%result)
         call Park_Field(field)
         if (allocated(result%r)) call move_alloc(result%r, field%r)
         if (allocated(result%z)) call move_alloc(result%z, field%z)
      end associate
   end subroutine End_Block_Transfer

end module ISFFT_Communication
