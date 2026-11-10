!------------------------------------------------------------------------------
! ISFFT_Allocation
! Preparation-time capacity reservation and final execution workspace setup.
!------------------------------------------------------------------------------
module ISFFT_Allocation
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_INT64_KIND
   use ISFFT_Plans, only: Axis_Plan, Kernel_Plan, poisson_plan, public_export, &
      public_import, public_plans, redistribution_forward_plan, redistribution_inverse_plan, &
      Transfer_Plan, Transform_Plan
   use ISFFT_Problem_State, only: dimensions, halo_width, internal_layout, ISFFT_PLAN_ESTIMATE, &
      Layout, original_layout, planning_level, public_fft_modes, &
      public_layout
   use ISFFT_Workspace, only: axis_imag, axis_real, axis_scratch, block_counts, &
      block_displs, block_recv_r, block_recv_z, block_send_r, &
      block_send_z, block_tile_r, block_tile_z, Clear_Field, &
      Destroy_Poisson_Operator, fft_precompute_slots, Field_Buffer, halo_work, &
      poisson_field, poisson_symbol, recv_complex_buffer, recv_real_buffer, &
      redistribution_field, Release_Execution_Scratch, Release_Layout_Workspace, Require_FFT_Precompute, &
      Reserve_Complex, Reserve_R, Reserve_Real, Reserve_Z, &
      Seal_Layout_Workspace, send_complex_buffer, send_real_buffer, Spectral_Field, &
      workspace_growths
   use ISFFT_MPI_Part, only: Layout_Box, Layout_Shape, Same_Layout
   use ISFFT_Communication, only: Release_Transfer, Release_Transform_Plan
   use ISFFT_Transform_Execution, only: Can_Direct_Input, Can_Direct_Output
   implicit none
   private
   public :: Finalize_Plans, Prepare_Execution_Workspaces, Reserve_Axis_Scratch, Reserve_Block, Reserve_FFT_Axis
   public :: Reserve_FFT_Base, Reserve_Field, Reserve_Halo_Workspace, Reserve_Pair, Reserve_R2C_Workspace
   public :: Reserve_Real_Shape, Reserve_Transfer, Reserve_Transfer_Workspace, Reserve_Transform

contains

   !------------------------------------------------------------------------------
   !! Problem preparation and lifecycle
   !------------------------------------------------------------------------------
   subroutine Prepare_Execution_Workspaces()
      integer :: which, mode, rank, ierr, n(3), np(3)
      logical :: direct

      !! Discard tuning buffers and reserve only the winning execution plans.
      call Clear_Field(poisson_field)
      call Clear_Field(redistribution_field)
      call Release_Execution_Scratch()
      call Release_Layout_Workspace()
      call Reserve_Transfer_Workspace(redistribution_forward_plan)
      call Reserve_Transfer_Workspace(redistribution_inverse_plan)
      call Reserve_Transfer(redistribution_field, redistribution_forward_plan)
      call Reserve_Transfer(redistribution_field, redistribution_inverse_plan)
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      do mode = 1, 2
         n = Layout_Shape(original_layout, rank)
         np = original_layout%process_grid
         if (mode == 2) then
            if (Same_Layout(original_layout, internal_layout)) cycle
            n = Layout_Shape(internal_layout, rank)
            np = internal_layout%process_grid
         end if
         if (halo_width > 0) then
            if (.not. Same_Layout(original_layout, internal_layout)) &
               call Reserve_Real_Shape(redistribution_field, n + 2 * halo_width, 1 - halo_width)
            call Reserve_Halo_Workspace(n, np, dimensions, halo_width)
         end if
      end do
      redistribution_field%prepared = .true.
      if (poisson_plan%ready) then
         direct = halo_width == 0 .and. Same_Layout(original_layout, internal_layout) .and. &
            Can_Direct_Input(poisson_plan) &
            .and. &
            Can_Direct_Output(poisson_plan)
         if (.not. direct) call Reserve_Field(poisson_field, original_layout)
         call Reserve_Transfer(poisson_field, redistribution_forward_plan)
         call Reserve_Transfer(poisson_field, redistribution_inverse_plan)
         call Reserve_Transform(poisson_field, poisson_plan, direct)
      end if
      do which = 1, 2
         if (.not. public_plans(which)%ready) cycle
         call Reserve_Transform(poisson_field, public_plans(which), Can_Direct_Input(public_plans(which)) .and. &
            Can_Direct_Output(public_plans(which)))
         do mode = 1, 2
            if (.not. btest(public_fft_modes(which), mode - 1)) cycle
            call Reserve_Field(poisson_field, public_layout(which, mode))
            call Reserve_Transfer(poisson_field, public_export(which, mode))
            call Reserve_Transfer(poisson_field, public_import(which, mode))
         end do
      end do
      poisson_field%prepared = .true.
      call Seal_Layout_Workspace()
   end subroutine ! Prepare_Execution_Workspaces

   !------------------------------------------------------------------------------
   !! Field and transform workspace preparation
   !------------------------------------------------------------------------------
   subroutine Reserve_Real_Shape(field, local_shape, lower)
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: local_shape(3), lower
      type(Field_Buffer), pointer :: node

      node => field%pool
      do while (associated(node))
         if (allocated(node%r)) then
            if (all(shape(node%r) == local_shape) .and. all(lbound(node%r) == lower)) return
         end if
         node => node%next
      end do
      allocate (node)
      allocate (node%r(lower:lower + local_shape(1) - 1, lower:lower + local_shape(2) - 1, lower:lower + &
         local_shape(3) - 1))
      node%next => field%pool
      field%pool => node
      workspace_growths = workspace_growths + 1
   end subroutine ! Reserve_Real_Shape

   subroutine Reserve_Field(field, field_layout, copies)
      type(Spectral_Field), intent(inout) :: field
      type(Layout), intent(in) :: field_layout
      integer, intent(in), optional :: copies
      type(Field_Buffer), pointer :: node
      integer :: rank, ierr, n(3), found, want

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      n = Layout_Shape(field_layout, rank)
      found = 0
      want = 1
      if (present(copies)) want = copies
      node => field%pool
      do while (associated(node))
         if (field_layout%complex_data) then
            if (allocated(node%z)) then
               if (all(shape(node%z) == n)) found = found + 1
            end if
         else
            if (allocated(node%r)) then
               if (all(shape(node%r) == n)) found = found + 1
            end if
         end if
         node => node%next
      end do
      do while (found < want)
         allocate (node)
         if (field_layout%complex_data) then
            allocate (node%z(n(1), n(2), n(3)))
         else
            allocate (node%r(n(1), n(2), n(3)))
         end if
         node%next => field%pool
         field%pool => node
         found = found + 1
         workspace_growths = workspace_growths + 1
      end do
   end subroutine ! Reserve_Field

   !! Equal shapes can still need two live buffers. Reserve both sides of a transfer.
   subroutine Reserve_Pair(field, source_layout, destination_layout)
      type(Spectral_Field), intent(inout) :: field
      type(Layout), intent(in) :: source_layout, destination_layout
      integer :: rank, ierr

      call Reserve_Field(field, source_layout)
      call Reserve_Field(field, destination_layout)
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      if (source_layout%complex_data .eqv. destination_layout%complex_data) then
         if (all(Layout_Shape(source_layout, rank) == Layout_Shape(destination_layout, rank))) &
            call Reserve_Field(field, source_layout, 2)
      end if
   end subroutine ! Reserve_Pair

   subroutine Reserve_Transfer(field, plan)
      type(Spectral_Field), intent(inout) :: field
      type(Transfer_Plan), intent(in) :: plan

      if (.not. plan%ready) return
      if (Same_Layout(plan%source, plan%destination)) return
      call Reserve_Pair(field, plan%source, plan%destination)
      call Reserve_Transfer_Workspace(plan)
   end subroutine ! Reserve_Transfer

   subroutine Reserve_Block(plan)
      type(Axis_Plan), intent(in) :: plan
      integer :: rank, ierr, n(3), nout(3), offset(3), parts
      integer(ISFFT_INT64_KIND) :: count, capacity, send_count
      logical :: need_real, need_complex

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call Layout_Box(plan%before%source, rank, n, offset)
      call Layout_Box(plan%after%destination, rank, nout, offset)
      capacity = max(1, plan%communication_batch_width)
      count = capacity * &
         max(plan%before%source%global_shape(plan%physical_axis), &
         plan%after%destination%global_shape(plan%physical_axis))
      send_count = max(1_ISFFT_INT64_KIND, capacity * max(n(plan%physical_axis), nout(plan%physical_axis)))
      need_real = .not. plan%before%source%complex_data .or. .not. plan%after%destination%complex_data
      need_complex = plan%before%source%complex_data .or. plan%after%destination%complex_data
      if (need_real) then
         call Reserve_Real(block_send_r, send_count)
         call Reserve_Real(block_recv_r, count)
         call Reserve_Real(block_tile_r, count)
      end if
      if (need_complex) then
         call Reserve_Complex(block_send_z, send_count)
         call Reserve_Complex(block_recv_z, count)
         call Reserve_Complex(block_tile_z, count)
      end if
      call MPI_Comm_size(ISFFT_Comm, parts, ierr)
      if (.not. allocated(block_counts)) allocate (block_counts(0:parts - 1), block_displs(0:parts - 1))
   end subroutine ! Reserve_Block

   subroutine Reserve_Transform(field, plan, direct_io)
      type(Spectral_Field), intent(inout) :: field
      type(Transform_Plan), intent(in) :: plan
      logical, intent(in), optional :: direct_io
      logical :: direct
      integer :: a

      if (.not. plan%ready) return
      direct = .false.
      if (present(direct_io)) direct = direct_io
      if (.not. direct) call Reserve_Field(field, plan%physical)
      call Reserve_Field(field, plan%spectral)
      do a = 1, plan%dimensions
         call Reserve_Axis_Scratch(plan%axes(a)%kernel)
         if (plan%axes(a)%tiled_block) then
            !! Tail batches have independently selected widths; reserve their scratch too.
            if (plan%axes(a)%tail_kernel%transform_length > 0) call Reserve_Axis_Scratch(plan%axes(a)%tail_kernel)
            call Reserve_Block(plan%axes(a))
            call Reserve_Pair(field, plan%axes(a)%before%source, plan%axes(a)%after%destination)
            cycle
         end if
         call Reserve_Transfer(field, plan%axes(a)%before)
         call Reserve_Transfer(field, plan%axes(a)%after)
         call Reserve_Transfer(field, plan%axes(a)%pack)
         call Reserve_Transfer(field, plan%axes(a)%unpack)
         call Reserve_Transfer(field, plan%axes(a)%reverse_before)
         call Reserve_Transfer(field, plan%axes(a)%reverse_after)
         call Reserve_Transfer(field, plan%axes(a)%reverse_pack)
         call Reserve_Transfer(field, plan%axes(a)%reverse_unpack)
         if (plan%axes(a)%transform_kind == 1) then
            if (a == 1 .and. direct) then
               call Reserve_Field(field, plan%axes(a)%output)
            else
               call Reserve_Pair(field, plan%axes(a)%input, plan%axes(a)%output)
            end if
         end if
      end do
   end subroutine ! Reserve_Transform

   !------------------------------------------------------------------------------
   !! Numerical workspace preparation
   !------------------------------------------------------------------------------
   subroutine Reserve_FFT_Base(slot)
      integer, intent(in) :: slot
      integer(ISFFT_INT64_KIND) :: n, m

      associate (p => fft_precompute_slots(slot))
         n = p%transform_length
         m = max(n, int(p%bluestein_m, ISFFT_INT64_KIND))
         call Reserve_Complex(p%line_work, n)
         call Reserve_Complex(p%split_work, 16 * n)
         call Reserve_Real(p%scalar_re, n)
         call Reserve_Real(p%scalar_im, n)
         if (.not. p%batch_ready) return
         call Reserve_Complex(p%batch_work, 16 * m)
         if (p%num_stages > 0) then
            call Reserve_Complex(p%mixed_1, 16 * n)
            call Reserve_Complex(p%mixed_2, 16 * n)
            call Reserve_Complex(p%mixed_3, 16 * n)
         end if
         if (p%rader_ready) call Reserve_Complex(p%rader_batch, 16 * max(1_ISFFT_INT64_KIND, n - 1))
      end associate
   end subroutine ! Reserve_FFT_Base

   subroutine Reserve_FFT_Axis(transform_length)
      integer, intent(in) :: transform_length
      integer :: slot

      call Require_FFT_Precompute(transform_length, slot)
      call Reserve_Complex(fft_precompute_slots(slot)%axis_batch, 16_ISFFT_INT64_KIND * transform_length)
      call Reserve_Complex(fft_precompute_slots(slot)%axis_lines, 16_ISFFT_INT64_KIND * transform_length)
   end subroutine ! Reserve_FFT_Axis

   subroutine Reserve_R2C_Workspace(transform_length, line_count, method)
      integer, intent(in) :: transform_length, line_count
      integer, intent(in), optional :: method
      integer :: slot, m, halfslot

      if (planning_level == ISFFT_PLAN_ESTIMATE .and. .not. present(method)) return
      if (present(method)) then
         if (method == 1) return
      end if
      call Require_FFT_Precompute(transform_length, slot)
      if (mod(transform_length, 2) == 0) then
         m = transform_length / 2
         call Reserve_Complex(fft_precompute_slots(slot)%r2c_pack, 16_ISFFT_INT64_KIND * m)
         call Require_FFT_Precompute(m, halfslot)
         if (.not. fft_precompute_slots(halfslot)%is_radix2) then
            call Reserve_Complex(fft_precompute_slots(slot)%r2c_natural, 16_ISFFT_INT64_KIND * m)
            call Reserve_Complex(fft_precompute_slots(slot)%r2c_generic, 16_ISFFT_INT64_KIND * m)
         end if
      end if
      if (mod(transform_length, 2) /= 0 .or. .not. fft_precompute_slots(slot)%is_radix2) &
         call Reserve_Complex(fft_precompute_slots(slot)%r2c_lines, int(transform_length, ISFFT_INT64_KIND) * &
         max(16, min(64, (line_count + 1) / 2)))
   end subroutine ! Reserve_R2C_Workspace

   subroutine Reserve_Axis_Scratch(plan)
      type(Kernel_Plan), intent(in) :: plan

      if (plan%transform_kind >= 3 .and. plan%complex_input) then
         call Reserve_Real(axis_real, int(plan%transform_length, ISFFT_INT64_KIND) * &
            max(1, min(plan%outer_batch_width, plan%line_count)))
         call Reserve_Real(axis_imag, int(plan%transform_length, ISFFT_INT64_KIND) * &
            max(1, min(plan%outer_batch_width, plan%line_count)))
      end if
      if (plan%memory_axis == 1) return
      if (plan%transform_kind == 2 .and. plan%direct_access) then
         call Reserve_FFT_Axis(plan%transform_length)
         return
      end if
      if (plan%transform_kind == 1 .or. (plan%transform_kind >= 3 .and. .not. plan%complex_input)) return
      call Reserve_Complex(axis_scratch, int(plan%transform_length, ISFFT_INT64_KIND) * &
         max(1, min(plan%outer_batch_width, plan%line_count)))
   end subroutine ! Reserve_Axis_Scratch

   !------------------------------------------------------------------------------
   !! Halo workspace preparation
   !------------------------------------------------------------------------------
   subroutine Reserve_Halo_Workspace(local_shape, process_grid, dimensions, halo_width)
      integer, intent(in) :: local_shape(3), process_grid(3), dimensions, halo_width
      integer :: a
      integer(ISFFT_INT64_KIND) :: capacity

      capacity = 0
      if (halo_width <= 0) return
      do a = 1, dimensions
         if (process_grid(a) <= 1) cycle
         capacity = max(capacity, int(halo_width, ISFFT_INT64_KIND) * product(int(local_shape, ISFFT_INT64_KIND)) / &
            max(1, local_shape(a)))
      end do
      if (capacity == 0) return
      if (allocated(halo_work)) then
         if (size(halo_work, 1, kind=ISFFT_INT64_KIND) >= capacity) return
         deallocate (halo_work)
      end if
      allocate (halo_work(capacity, 4))
      workspace_growths = workspace_growths + 1
   end subroutine ! Reserve_Halo_Workspace

   subroutine Reserve_Transfer_Workspace(plan)
      type(Transfer_Plan), intent(in) :: plan
      integer(ISFFT_INT64_KIND) :: send_count, receive_count

      if (.not. plan%ready) return
      if (Same_Layout(plan%source, plan%destination) .or. plan%local_only) return
      if ((plan%direct .or. plan%sparse) .and. .not. (plan%chunked .or. plan%streamed) .and. &
         plan%gather_comm == MPI_COMM_NULL) return
      if (plan%chunked .or. plan%streamed) then
         send_count = max(1_ISFFT_INT64_KIND, min(plan%chunk_limit, plan%max_peer, plan%stream_width))
         receive_count = send_count
      else
         send_count = max(1_ISFFT_INT64_KIND, sum(plan%send_counts))
         receive_count = max(1_ISFFT_INT64_KIND, sum(plan%recv_counts))
         if (plan%gather_comm /= MPI_COMM_NULL) send_count = max(1_ISFFT_INT64_KIND, plan%gather_count)
      end if
      if (plan%source%complex_data) then
         call Reserve_Z(send_complex_buffer, send_count)
         call Reserve_Z(recv_complex_buffer, receive_count)
      else
         call Reserve_R(send_real_buffer, send_count)
         call Reserve_R(recv_real_buffer, receive_count)
      end if
   end subroutine ! Reserve_Transfer_Workspace

   subroutine Finalize_Plans()
      integer :: i, m

      call Clear_Field(poisson_field)
      call Clear_Field(redistribution_field)
      call Release_Transfer(redistribution_forward_plan)
      call Release_Transfer(redistribution_inverse_plan)
      call Release_Transform_Plan(poisson_plan)
      do i = 1, 2
         call Release_Transform_Plan(public_plans(i))
         do m = 1, 2
            call Release_Transfer(public_export(i, m))
            call Release_Transfer(public_import(i, m))
         end do
      end do
      call Release_Layout_Workspace()
      call Destroy_Poisson_Operator(poisson_symbol)
   end subroutine ! Finalize_Plans

end module ISFFT_Allocation
