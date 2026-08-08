module ISFFT_2D_Mixed_Poisson
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_FFT_kernel, only: ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c
   use ISFFT_DST_kernel, only: ISFFT_1D_DST1
   use ISFFT_DCT_kernel, only: ISFFT_1D_DCT1
   use ISFFT_Discrete_Symbols
   use ISFFT_LocalTranspose, only: local_transpose_01_real, local_transpose_01_complex
   use ISFFT_GlobalTranspose, only: global_transpose_01_real, global_transpose_01_complex, &
      global_transpose_10_real, global_transpose_10_complex
   implicit none

   private
   public :: ISFFT_2D_Mixed_Poisson_serial, ISFFT_2D_Mixed_Poisson_block
   public :: ISFFT_2D_Mixed_Poisson_slab

contains

   subroutine ISFFT_2D_Mixed_Poisson_serial(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: &
         phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(in) :: &
         f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND) :: work_re(nx_global,ny_global,1)
      real(kind=ISFFT_REAL_KIND) :: work_im(nx_global,ny_global,1)

      work_re = f(1:nx_global,1:ny_global,1:1)
      work_im = 0.0_ISFFT_REAL_KIND
      call Transform_X_Physical(work_re, work_im, +1)
      call Transform_Y_Physical(work_re, work_im, +1)
      call Algebra_XY(work_re, work_im, 0, 0)
      call Transform_Y_Physical(work_re, work_im, -1)
      call Transform_X_Physical(work_re, work_im, -1)

      phi = 0.0_ISFFT_REAL_KIND
      phi(1:nx_global,1:ny_global,1) = work_re(:,:,1)
      call Apply_Dirichlet_Boundaries(phi)
   end subroutine ISFFT_2D_Mixed_Poisson_serial

   subroutine ISFFT_2D_Mixed_Poisson_block(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: &
         phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(in) :: &
         f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re(1:nx,1:ny,1) = f(1:nx,1:ny,1)

      call Transform_Axis_Block(1, +1, f_re, f_im)
      call Transform_Axis_Block(2, +1, f_re, f_im)
      call Algebra_XY(f_re, f_im, i_offset(npx)-1, j_offset(npy)-1)
      call Transform_Axis_Block(2, -1, f_re, f_im)
      call Transform_Axis_Block(1, -1, f_re, f_im)

      phi = 0.0_ISFFT_REAL_KIND
      phi(1:nx,1:ny,1) = f_re(1:nx,1:ny,1)
      call Apply_Dirichlet_Boundaries(phi)
   end subroutine ISFFT_2D_Mixed_Poisson_block

   subroutine ISFFT_2D_Mixed_Poisson_slab(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: &
         phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(in) :: &
         f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      integer :: x_offset_layout
      logical :: complex_after_x

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re(1:nx,1:ny,1) = f(1:nx,1:ny,1)

      call Transform_X_Physical(f_re, f_im, +1)
      complex_after_x = (BC(1) == ISFFT_BC_PERIODIC)
      if (complex_after_x) then
         call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
         call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
      else
         call local_transpose_01_real(f_re, phi_re_transposed_a)
         call global_transpose_01_real(phi_re_transposed_a, f_hat_y)
         f_hat_im_y = 0.0_ISFFT_REAL_KIND
      end if

      call Transform_Y_Layout(f_hat_y, f_hat_im_y, +1)
      x_offset_layout = npy*size(f_hat_y,2)
      call Algebra_YX(f_hat_y, f_hat_im_y, x_offset_layout)
      call Transform_Y_Layout(f_hat_y, f_hat_im_y, -1)

      if (complex_after_x) then
         call global_transpose_10_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
         call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, phi_re, phi_im)
      else
         call global_transpose_10_real(f_hat_y, phi_re_transposed_a)
         call local_transpose_01_real(phi_re_transposed_a, phi_re)
         phi_im = 0.0_ISFFT_REAL_KIND
      end if
      call Transform_X_Physical(phi_re, phi_im, -1)

      phi = 0.0_ISFFT_REAL_KIND
      phi(1:nx,1:ny,1) = phi_re(1:nx,1:ny,1)
      call Apply_Dirichlet_Boundaries(phi)
   end subroutine ISFFT_2D_Mixed_Poisson_slab

   subroutine Transform_X_Physical(a_re, a_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      integer :: j
      logical :: transform_im

      transform_im = Input_Imag_Active(1, sign)
      do j = 1, size(a_re,2)
         call Transform_Line(BC(1), 1, nx_global, a_re(:,j,1), a_im(:,j,1), sign, transform_im)
      end do
   end subroutine Transform_X_Physical

   subroutine Transform_Y_Physical(a_re, a_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      integer :: i, i_stop
      logical :: transform_im

      i_stop = size(a_re,1)
      if (First_Periodic_Axis() == 1) i_stop = min(i_stop, halfspectrum_x)
      transform_im = Input_Imag_Active(2, sign)
      do i = 1, i_stop
         call Transform_Line(BC(2), 2, ny_global, a_re(i,:,1), a_im(i,:,1), sign, transform_im)
      end do
   end subroutine Transform_Y_Physical

   subroutine Transform_Y_Layout(a_re, a_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      integer :: i, i_stop, x_offset
      logical :: transform_im

      x_offset = npy*size(a_re,2)
      i_stop = size(a_re,2)
      if (First_Periodic_Axis() == 1) &
         i_stop = max(0, min(i_stop, halfspectrum_x-x_offset))
      transform_im = Input_Imag_Active(2, sign)
      do i = 1, i_stop
         call Transform_Line(BC(2), 2, ny_global, a_re(:,i,1), a_im(:,i,1), sign, transform_im)
      end do
   end subroutine Transform_Y_Layout

   subroutine Transform_Axis_Block(axis, sign, a_re, a_im)
      integer, intent(in) :: axis, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer :: i, j, ierr, gi, gj, i_stop, j_stop
      logical :: transform_im, output_im

      transform_im = Input_Imag_Active(axis, sign)
      output_im = Output_Imag_Active(axis, sign)
      if (axis == 1) then
         j_stop = Active_Local_Stop(2, axis)
         do j = 1, j_stop
            x_sendbuf = 0.0_ISFFT_REAL_KIND
            x_im_sendbuf = 0.0_ISFFT_REAL_KIND
            do i = 1, nx
               gi = i_offset(npx)+i-1
               x_sendbuf(gi) = a_re(i,j,1)
               if (transform_im) x_im_sendbuf(gi) = a_im(i,j,1)
            end do
            call MPI_Allreduce(x_sendbuf, x_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
            if (transform_im) then
               call MPI_Allreduce(x_im_sendbuf, x_im_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
            else
               x_im_recvbuf = 0.0_ISFFT_REAL_KIND
            end if
            call Transform_Line(BC(1), 1, nx_global, x_recvbuf, x_im_recvbuf, sign, transform_im)
            do i = 1, nx
               gi = i_offset(npx)+i-1
               a_re(i,j,1) = x_recvbuf(gi)
               if (output_im) a_im(i,j,1) = x_im_recvbuf(gi)
            end do
         end do
      else
         i_stop = Active_Local_Stop(1, axis)
         do i = 1, i_stop
            y_sendbuf = 0.0_ISFFT_REAL_KIND
            y_im_sendbuf = 0.0_ISFFT_REAL_KIND
            do j = 1, ny
               gj = j_offset(npy)+j-1
               y_sendbuf(gj) = a_re(i,j,1)
               if (transform_im) y_im_sendbuf(gj) = a_im(i,j,1)
            end do
            call MPI_Allreduce(y_sendbuf, y_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
            if (transform_im) then
               call MPI_Allreduce(y_im_sendbuf, y_im_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
            else
               y_im_recvbuf = 0.0_ISFFT_REAL_KIND
            end if
            call Transform_Line(BC(2), 2, ny_global, y_recvbuf, y_im_recvbuf, sign, transform_im)
            do j = 1, ny
               gj = j_offset(npy)+j-1
               a_re(i,j,1) = y_recvbuf(gj)
               if (output_im) a_im(i,j,1) = y_im_recvbuf(gj)
            end do
         end do
      end if
   end subroutine Transform_Axis_Block

   subroutine Transform_Line(bc, axis, n, line_re, line_im, sign, transform_im)
      integer, intent(in) :: bc, axis, n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(n), line_im(n)
      logical, intent(in) :: transform_im

      select case (bc)
      case (ISFFT_BC_PERIODIC)
         if (axis == First_Periodic_Axis()) then
            call ISFFT_1D_FFT_r2c(n, line_re, line_im, sign)
         else
            call ISFFT_1D_FFT_c2c(n, line_re, line_im, sign)
         end if
      case (ISFFT_BC_DIRICHLET)
         if (n > 2) then
            call ISFFT_1D_DST1(n-2, line_re(2:n-1), sign)
            if (transform_im) call ISFFT_1D_DST1(n-2, line_im(2:n-1), sign)
         end if
         line_re(1) = 0.0_ISFFT_REAL_KIND
         line_re(n) = 0.0_ISFFT_REAL_KIND
         if (transform_im) then
            line_im(1) = 0.0_ISFFT_REAL_KIND
            line_im(n) = 0.0_ISFFT_REAL_KIND
         end if
      case (ISFFT_BC_NEUMANN)
         call ISFFT_1D_DCT1(n, line_re, sign)
         if (transform_im) call ISFFT_1D_DCT1(n, line_im, sign)
      end select
   end subroutine Transform_Line

   subroutine Algebra_XY(a_re, a_im, x_offset, y_offset)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: x_offset, y_offset
      integer :: i, j
      real(kind=ISFFT_REAL_KIND), allocatable :: lambda_x(:), lambda_y(:)
      logical, allocatable :: boundary_x(:), boundary_y(:), active_x(:), active_y(:)

      allocate(lambda_x(size(a_re,1)), lambda_y(size(a_re,2)))
      allocate(boundary_x(size(a_re,1)), boundary_y(size(a_re,2)))
      allocate(active_x(size(a_re,1)), active_y(size(a_re,2)))
      call Precompute_Modes(1, size(a_re,1), x_offset, lambda_x, boundary_x, active_x)
      call Precompute_Modes(2, size(a_re,2), y_offset, lambda_y, boundary_y, active_y)

      do j = 1, size(a_re,2)
         do i = 1, size(a_re,1)
            call Divide_Mode(a_re(i,j,1), a_im(i,j,1), lambda_x(i)+lambda_y(j), &
               active_x(i) .and. active_y(j) .and. .not. boundary_x(i) .and. .not. boundary_y(j))
         end do
      end do
      deallocate(lambda_x, lambda_y, boundary_x, boundary_y, active_x, active_y)
   end subroutine Algebra_XY

   subroutine Algebra_YX(a_re, a_im, x_offset)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: x_offset
      integer :: i, j
      real(kind=ISFFT_REAL_KIND), allocatable :: lambda_x(:), lambda_y(:)
      logical, allocatable :: boundary_x(:), boundary_y(:), active_x(:), active_y(:)

      allocate(lambda_x(size(a_re,2)), lambda_y(size(a_re,1)))
      allocate(boundary_x(size(a_re,2)), boundary_y(size(a_re,1)))
      allocate(active_x(size(a_re,2)), active_y(size(a_re,1)))
      call Precompute_Modes(1, size(a_re,2), x_offset, lambda_x, boundary_x, active_x)
      call Precompute_Modes(2, size(a_re,1), 0, lambda_y, boundary_y, active_y)

      do i = 1, size(a_re,2)
         do j = 1, size(a_re,1)
            call Divide_Mode(a_re(j,i,1), a_im(j,i,1), lambda_x(i)+lambda_y(j), &
               active_x(i) .and. active_y(j) .and. .not. boundary_x(i) .and. .not. boundary_y(j))
         end do
      end do
      deallocate(lambda_x, lambda_y, boundary_x, boundary_y, active_x, active_y)
   end subroutine Algebra_YX

   subroutine Precompute_Modes(axis, local_n, offset, lambda, boundary, active)
      integer, intent(in) :: axis, local_n, offset
      real(kind=ISFFT_REAL_KIND), intent(out) :: lambda(local_n)
      logical, intent(out) :: boundary(local_n), active(local_n)
      type(ISFFT_Scheme_Coefficients) :: coeff
      integer :: index, mode, n_global
      real(kind=ISFFT_REAL_KIND) :: length

      call ISFFT_Get_Scheme_Coefficients(If_scheme, coeff)
      if (axis == 1) then
         n_global = nx_global
         length = SLx
      else
         n_global = ny_global
         length = SLy
      end if
      do index = 1, local_n
         call Mode_Info(axis, index+offset, mode, boundary(index), active(index))
         lambda(index) = ISFFT_Lambda_For_BC(BC(axis), mode, n_global, length, If_scheme, coeff)
      end do
   end subroutine Precompute_Modes

   subroutine Divide_Mode(value_re, value_im, lambda, active)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: value_re, value_im
      real(kind=ISFFT_REAL_KIND), intent(in) :: lambda
      logical, intent(in) :: active

      if (active .and. lambda /= 0.0_ISFFT_REAL_KIND) then
         value_re = -value_re/lambda
         value_im = -value_im/lambda
      else
         value_re = 0.0_ISFFT_REAL_KIND
         value_im = 0.0_ISFFT_REAL_KIND
      end if
   end subroutine Divide_Mode

   subroutine Mode_Info(axis, global_index, mode, on_boundary, active)
      integer, intent(in) :: axis, global_index
      integer, intent(out) :: mode
      logical, intent(out) :: on_boundary, active
      integer :: n_global, half_limit

      n_global = merge(nx_global, ny_global, axis == 1)
      half_limit = n_global/2+1
      active = (axis /= First_Periodic_Axis() .or. global_index <= half_limit)
      on_boundary = .false.
      mode = global_index-1
      select case (BC(axis))
      case (ISFFT_BC_PERIODIC)
         if (2*mode >= n_global) mode = mode-n_global
      case (ISFFT_BC_DIRICHLET)
         on_boundary = (global_index == 1 .or. global_index == n_global)
         if (on_boundary) mode = 0
      end select
   end subroutine Mode_Info

   integer function First_Periodic_Axis() result(axis)
      axis = 0
      if (BC(1) == ISFFT_BC_PERIODIC) then
         axis = 1
      else if (BC(2) == ISFFT_BC_PERIODIC) then
         axis = 2
      end if
   end function First_Periodic_Axis

   logical function Input_Imag_Active(axis, sign)
      integer, intent(in) :: axis, sign
      if (sign > 0) then
         Input_Imag_Active = any(BC(1:axis-1) == ISFFT_BC_PERIODIC)
      else
         Input_Imag_Active = any(BC(1:axis) == ISFFT_BC_PERIODIC)
      end if
   end function Input_Imag_Active

   logical function Output_Imag_Active(axis, sign)
      integer, intent(in) :: axis, sign
      if (sign > 0) then
         Output_Imag_Active = any(BC(1:axis) == ISFFT_BC_PERIODIC)
      else if (axis > 1) then
         Output_Imag_Active = any(BC(1:axis-1) == ISFFT_BC_PERIODIC)
      else
         Output_Imag_Active = .false.
      end if
   end function Output_Imag_Active

   integer function Active_Local_Stop(dim, transform_axis) result(stop_idx)
      integer, intent(in) :: dim, transform_axis
      integer :: local_n, local_offset, half_limit

      if (dim == 1) then
         local_n = nx
         local_offset = i_offset(npx)
         half_limit = halfspectrum_x
      else
         local_n = ny
         local_offset = j_offset(npy)
         half_limit = halfspectrum_y
      end if
      if (First_Periodic_Axis() == dim .and. transform_axis > dim) then
         stop_idx = max(0, min(local_n, half_limit-local_offset+1))
      else
         stop_idx = local_n
      end if
   end function Active_Local_Stop

   subroutine Apply_Dirichlet_Boundaries(phi)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: &
         phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      integer :: i, j, gi, gj

      if (BC(1) == ISFFT_BC_DIRICHLET) then
         do j = 1, ny
            do i = 1, nx
               gi = i_offset(npx)+i-1
               if (gi == 1 .or. gi == nx_global) phi(i,j,1) = 0.0_ISFFT_REAL_KIND
            end do
         end do
      end if
      if (BC(2) == ISFFT_BC_DIRICHLET) then
         do j = 1, ny
            gj = j_offset(npy)+j-1
            if (gj == 1 .or. gj == ny_global) phi(1:nx,j,1) = 0.0_ISFFT_REAL_KIND
         end do
      end if
   end subroutine Apply_Dirichlet_Boundaries

end module ISFFT_2D_Mixed_Poisson
