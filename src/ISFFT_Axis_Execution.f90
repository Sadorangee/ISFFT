!------------------------------------------------------------------------------
! ISFFT_Axis_Execution
! Local single-axis access, batching and numerical-kernel dispatch.
!------------------------------------------------------------------------------
module ISFFT_Axis_Execution
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND
   use ISFFT_Plans, only: Axis_Plan, Kernel_Plan
   use ISFFT_Workspace, only: axis_imag, axis_real, axis_scratch, Get_Complex, &
      Get_Real, Park_Field, Spectral_Field
   use ISFFT_MPI_Part, only: Axis_Position, Layout_Shape
   use ISFFT_DCT_kernel, only: ISFFT_Execute_DCT_Plan
   use ISFFT_DST_kernel, only: ISFFT_Execute_DST_Plan
   use ISFFT_FFT_kernel, only: ISFFT_Execute_Axis_Lines, ISFFT_Execute_Complex_Plan, &
      ISFFT_Execute_R2C_Forward, ISFFT_Execute_R2C_Plan
   implicit none
   private
   public :: Axis_Execute, Execute_Kernel, Execute_Kernel_Forward

   interface Execute_Kernel
      module procedure Execute_Kernel_R2C, Execute_Kernel_Real, Execute_Kernel_Complex
   end interface

contains

   !------------------------------------------------------------------------------
   !! Single-axis execution
   !------------------------------------------------------------------------------
   subroutine Axis_Execute(plan, field, sign)
      type(Axis_Plan), intent(in) :: plan
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: sign
      integer :: n(3), rank, ierr
      real(kind=ISFFT_REAL_KIND), allocatable :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: z(:, :, :)

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      if (plan%transform_kind == 1) then
         if (sign > 0) then
            n = Layout_Shape(plan%output, rank)
            call Get_Complex(field, n, z)
            if (size(field%r, kind=ISFFT_INT64_KIND) > 0) call Execute_Kernel(plan%kernel, field%r, z, 1)
            call Park_Field(field)
            call move_alloc(z, field%z)
         else
            n = Layout_Shape(plan%input, rank)
            call Get_Real(field, n, r)
            if (size(r, kind=ISFFT_INT64_KIND) > 0) call Execute_Kernel(plan%kernel, r, field%z, -1)
            call Park_Field(field)
            call move_alloc(r, field%r)
         end if
      else if (allocated(field%z)) then
         if (size(field%z, kind=ISFFT_INT64_KIND) > 0) &
            call Execute_Kernel(plan%kernel, field%z, sign, Axis_Position(plan%input%axis_order, plan%physical_axis))
      else
         if (size(field%r, kind=ISFFT_INT64_KIND) > 0) call Execute_Kernel(plan%kernel, field%r, sign)
      end if
   end subroutine ! Axis_Execute

   !------------------------------------------------------------------------------
   !! Outer grouping and numerical kernel dispatch
   !------------------------------------------------------------------------------
   !! Outer grouping width and kernel width are intentionally independent.
   subroutine Execute_Kernel_Forward(plan, re, spectrum)
      type(Kernel_Plan), intent(in) :: plan
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: z(:, :, :)
      integer :: lines, first, last, capacity

      lines = size(re, kind=ISFFT_INT64_KIND) / plan%transform_length
      if (lines == 0) return
      r(1:plan%transform_length, 1:lines, 1:1) => re
      z(1:plan%transform_length / 2 + 1, 1:lines, 1:1) => spectrum
      capacity = max(1, plan%outer_batch_width)
      do first = 1, lines, capacity
         last = min(lines, first + capacity - 1)
         call ISFFT_Execute_R2C_Forward(plan%transform_length, r(:, first:last, :), z(:, first:last, :), &
            plan%kernel_method, &
            plan%inner_batch_width)
      end do
   end subroutine ! Execute_Kernel_Forward

   subroutine Execute_Kernel_R2C(plan, re, spectrum, sign)
      type(Kernel_Plan), intent(in) :: plan
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: z(:, :, :)
      integer :: lines, first, last, capacity

      lines = size(re, kind=ISFFT_INT64_KIND) / plan%transform_length
      if (lines == 0) return
      r(1:plan%transform_length, 1:lines, 1:1) => re
      z(1:plan%transform_length / 2 + 1, 1:lines, 1:1) => spectrum
      capacity = max(1, plan%outer_batch_width)
      do first = 1, lines, capacity
         last = min(lines, first + capacity - 1)
         call ISFFT_Execute_R2C_Plan(plan%transform_length, r(:, first:last, :), z(:, first:last, :), sign, &
            plan%kernel_method, &
            plan%inner_batch_width)
      end do
   end subroutine ! Execute_Kernel_R2C

   subroutine Execute_Kernel_Real(plan, values, sign)
      type(Kernel_Plan), intent(in) :: plan
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: values(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: r(:, :, :)
      integer :: lines, first, last, capacity

      lines = size(values, kind=ISFFT_INT64_KIND) / plan%transform_length
      if (lines == 0) return
      r(1:plan%transform_length, 1:lines, 1:1) => values
      capacity = max(1, plan%outer_batch_width)
      do first = 1, lines, capacity
         last = min(lines, first + capacity - 1)
         if (plan%transform_kind == 3) then
            call ISFFT_Execute_DST_Plan(plan%transform_length, r(:, first:last, :), sign, plan%kernel_method, &
               plan%inner_batch_width)
         else
            call ISFFT_Execute_DCT_Plan(plan%transform_length, r(:, first:last, :), sign, plan%kernel_method, &
               plan%inner_batch_width)
         end if
      end do
   end subroutine ! Execute_Kernel_Real

   subroutine Execute_Complex_Group(plan, values, sign)
      type(Kernel_Plan), intent(in) :: plan
      integer, intent(in) :: sign
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous :: values(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: r(:, :, :), im(:, :, :)
      integer :: ns

      if (plan%transform_kind == 2) then
         call ISFFT_Execute_Complex_Plan(plan%transform_length, values, sign, plan%kernel_method, &
            plan%inner_batch_width)
      else
         ns = size(values)
         if (.not. allocated(axis_real)) error stop 'complex R2R workspace not prepared'
         if (size(axis_real) < ns) error stop 'complex R2R workspace too small'
         r(1:size(values, 1), 1:size(values, 2), 1:size(values, 3)) => axis_real(:ns)
         im(1:size(values, 1), 1:size(values, 2), 1:size(values, 3)) => axis_imag(:ns)
         r = real(values, ISFFT_REAL_KIND)
         im = aimag(values)
         if (plan%transform_kind == 3) then
            call ISFFT_Execute_DST_Plan(plan%transform_length, r, sign, plan%kernel_method, plan%inner_batch_width)
            call ISFFT_Execute_DST_Plan(plan%transform_length, im, sign, plan%kernel_method, plan%inner_batch_width)
         else
            call ISFFT_Execute_DCT_Plan(plan%transform_length, r, sign, plan%kernel_method, plan%inner_batch_width)
            call ISFFT_Execute_DCT_Plan(plan%transform_length, im, sign, plan%kernel_method, plan%inner_batch_width)
         end if
         values = cmplx(r, im, ISFFT_REAL_KIND)
      end if
   end subroutine ! Execute_Complex_Group

   subroutine Execute_Kernel_Complex(plan, values, sign, memory_axis, selected_lines)
      type(Kernel_Plan), intent(in) :: plan
      integer, intent(in) :: sign
      integer, intent(in), optional :: memory_axis
      integer, intent(in), optional :: selected_lines(:)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: values(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: z(:), group(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: packed(:, :, :)
      integer :: lines, first, last, capacity, a, line, id, i, j, k, width

      if (size(values, kind=ISFFT_INT64_KIND) == 0) return
      a = 1
      if (present(memory_axis)) a = memory_axis
      lines = size(values, kind=ISFFT_INT64_KIND) / plan%transform_length
      if (present(selected_lines)) lines = size(selected_lines)
      capacity = max(1, min(plan%outer_batch_width, lines))
      if (a == 1) then
         z(1:size(values)) => values
         do first = 1, lines, capacity
            last = min(lines, first + capacity - 1)
            group(1:plan%transform_length, 1:last - first + 1, 1:1) &
               => z((first - 1) * plan%transform_length + 1:last * plan%transform_length)
            call Execute_Complex_Group(plan, group, sign)
         end do
         return
      end if
      !! The existing direct kernel consumes the user's strided layout, fusing its
      !! loads/stores with FFT work. Selection happens only during initialization.
      if (plan%transform_kind == 2 .and. plan%direct_access) then
         call ISFFT_Execute_Axis_Lines(plan%transform_length, values, a, sign, plan%direct_batch_width, &
            planned_method=plan%kernel_method, &
            selected_lines=selected_lines)
         return
      end if
      if (.not. allocated(axis_scratch)) error stop 'axis workspace not prepared'
      if (size(axis_scratch, kind=ISFFT_INT64_KIND) < int(plan%transform_length, ISFFT_INT64_KIND) * capacity) &
         error stop 'axis workspace too small'
      packed(1:plan%transform_length, 1:capacity, 1:1) => axis_scratch(1:plan%transform_length * capacity)
      do first = 0, lines - 1, capacity
         width = min(capacity, lines - first)
         do line = 1, width
            id = first + line - 1
            if (present(selected_lines)) id = selected_lines(id + 1) - 1
            if (a == 2) then
               i = mod(id, size(values, 1)) + 1
               k = id / size(values, 1) + 1
               packed(:, line, 1) = values(i, :, k)
            else
               i = mod(id, size(values, 1)) + 1
               j = id / size(values, 1) + 1
               packed(:, line, 1) = values(i, j, :)
            end if
         end do
         group(1:plan%transform_length, 1:width, 1:1) => axis_scratch(1:plan%transform_length * width)
         call Execute_Complex_Group(plan, group, sign)
         do line = 1, width
            id = first + line - 1
            if (present(selected_lines)) id = selected_lines(id + 1) - 1
            if (a == 2) then
               i = mod(id, size(values, 1)) + 1
               k = id / size(values, 1) + 1
               values(i, :, k) = packed(:, line, 1)
            else
               i = mod(id, size(values, 1)) + 1
               j = id / size(values, 1) + 1
               values(i, j, :) = packed(:, line, 1)
            end if
         end do
      end do
   end subroutine ! Execute_Kernel_Complex

end module ISFFT_Axis_Execution
