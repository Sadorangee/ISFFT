module ISFFT_3D_Mixed_Poisson
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_FFT_kernel
   use ISFFT_DST_kernel
   use ISFFT_DCT_kernel
   use ISFFT_Discrete_Symbols
   use ISFFT_LocalTranspose
   use ISFFT_GlobalTranspose
   implicit none

   private
   public :: ISFFT_3D_Mixed_Poisson_serial, ISFFT_3D_Mixed_Poisson_block
   public :: ISFFT_3D_Mixed_Poisson_transpose

   integer, parameter :: MIXED_TRANSFORM_FFT_R2C = 1
   integer, parameter :: MIXED_TRANSFORM_FFT_C2C = 2
   integer, parameter :: MIXED_TRANSFORM_DST = 3
   integer, parameter :: MIXED_TRANSFORM_DCT = 4

contains

   subroutine ISFFT_3D_Mixed_Poisson_serial(phi, f)
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: f
      real(kind=ISFFT_REAL_KIND), dimension(1:nx,1:ny,1:nz) :: work_re, work_im
      integer :: i, j, k

      work_re = 0.0_ISFFT_REAL_KIND
      work_im = 0.0_ISFFT_REAL_KIND
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               work_re(i,j,k) = f(i,j,k)
            end do
         end do
      end do

      call Mixed_Forward_Serial(work_re, work_im)
      call Mixed_Algebra_Block(work_re, work_im)
      call Mixed_Inverse_Serial(work_re, work_im)

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               phi(i,j,k) = work_re(i,j,k)
            end do
         end do
      end do
      call Mixed_Dirichlet_HomoBC(phi)
   end subroutine ISFFT_3D_Mixed_Poisson_serial

   subroutine ISFFT_3D_Mixed_Poisson_block(phi, f)
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: f
      integer :: i, j, k, ierr

      f_re = 0.0_ISFFT_REAL_KIND
      if (Mixed_Any_Periodic()) f_im = 0.0_ISFFT_REAL_KIND
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               f_re(i,j,k) = f(i,j,k)
            end do
         end do
      end do

      call Mixed_Forward_Block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Mixed_Algebra_Block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Mixed_Inverse_Block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               phi(i,j,k) = f_re(i,j,k)
            end do
         end do
      end do
      call Mixed_Dirichlet_HomoBC(phi)
   end subroutine ISFFT_3D_Mixed_Poisson_block

   subroutine ISFFT_3D_Mixed_Poisson_transpose(phi, f)
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), intent(inout) :: f
      integer :: i, j, k, ierr

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               f_re(i,j,k) = f(i,j,k)
            end do
         end do
      end do

      call Mixed_Transform_X_Physical(f_re, f_im, +1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (Mixed_Output_Imag_Active(1, +1)) then
         if (If_decomp == 1) then
            call local_transpose_01_complex(f_re, f_im, f_hat_y, f_hat_im_y)
         else
            call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
            call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         end if
      else
         if (If_decomp == 1) then
            call local_transpose_01_real(f_re, f_hat_y)
         else
            call local_transpose_01_real(f_re, phi_re_transposed_a)
            call global_transpose_01_real(phi_re_transposed_a, f_hat_y)
         end if
      end if
      call Mixed_Transform_Y_Layout(f_hat_y, f_hat_im_y, +1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (Mixed_Output_Imag_Active(2, +1)) then
         if (If_decomp == 1) then
            call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
            call global_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat, f_hat_im)
         else
            call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_b, phi_im_transposed_b)
            call global_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, f_hat, f_hat_im)
         end if
      else
         if (If_decomp == 1) then
            call local_transpose_02_real(f_hat_y, phi_re_transposed_a)
            call global_transpose_02_real(phi_re_transposed_a, f_hat)
         else
            call local_transpose_02_real(f_hat_y, phi_re_transposed_b)
            call global_transpose_02_real(phi_re_transposed_b, f_hat)
         end if
      end if
      call Mixed_Transform_Z_Layout(f_hat, f_hat_im, +1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Mixed_Algebra_ZXY(f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Mixed_Transform_Z_Layout(f_hat, f_hat_im, -1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (Mixed_Output_Imag_Active(3, -1)) then
         if (If_decomp == 1) then
            call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_a, phi_im_transposed_a)
            call local_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, phi_hat_x, phi_hat_im_x)
         else
            call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_b, phi_im_transposed_b)
            call local_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, phi_hat_x, phi_hat_im_x)
         end if
      else
         if (If_decomp == 1) then
            call global_transpose_20_real(f_hat, phi_re_transposed_a)
            call local_transpose_02_real(phi_re_transposed_a, phi_hat_x)
         else
            call global_transpose_20_real(f_hat, phi_re_transposed_b)
            call local_transpose_02_real(phi_re_transposed_b, phi_hat_x)
         end if
      end if
      call Mixed_Transform_Y_Layout(phi_hat_x, phi_hat_im_x, -1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (Mixed_Output_Imag_Active(2, -1)) then
         if (If_decomp == 1) then
            call local_transpose_01_complex(phi_hat_x, phi_hat_im_x, phi_re, phi_im)
         else
            call global_transpose_10_complex(phi_hat_x, phi_hat_im_x, phi_re_transposed_a, phi_im_transposed_a)
            call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, phi_re, phi_im)
         end if
      else
         if (If_decomp == 1) then
            call local_transpose_01_real(phi_hat_x, phi_re)
         else
            call global_transpose_10_real(phi_hat_x, phi_re_transposed_a)
            call local_transpose_01_real(phi_re_transposed_a, phi_re)
         end if
      end if
      call Mixed_Transform_X_Physical(phi_re, phi_im, -1)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               phi(i,j,k) = phi_re(i,j,k)
            end do
         end do
      end do
      call Mixed_Dirichlet_HomoBC(phi)
   end subroutine ISFFT_3D_Mixed_Poisson_transpose

   subroutine Mixed_Forward_Block(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)

      call Mixed_Transform_Axis_Block(1, BC(1), +1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Block(2, BC(2), +1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Block(3, BC(3), +1, phi_re_inout, phi_im_inout)
   end subroutine Mixed_Forward_Block

   subroutine Mixed_Inverse_Block(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)

      call Mixed_Transform_Axis_Block(3, BC(3), -1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Block(2, BC(2), -1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Block(1, BC(1), -1, phi_re_inout, phi_im_inout)
   end subroutine Mixed_Inverse_Block

   subroutine Mixed_Transform_X_Physical(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      integer, intent(in) :: sign
      integer :: j, k
      logical :: transform_im

      transform_im = Mixed_Input_Imag_Active(1, sign)
      select case (BC(1))
      case (ISFFT_BC_PERIODIC)
         if (1 == Mixed_First_Periodic_Axis()) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call ISFFT_1D_FFT_r2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call ISFFT_1D_FFT_c2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_DIRICHLET)
         if (transform_im) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call Mixed_Line_DST_Complex(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call Mixed_Line_DST_Real(nx_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_NEUMANN)
         if (transform_im) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call Mixed_Line_DCT_Complex(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, size(phi_re_inout, 2)
                  call ISFFT_1D_DCT1(nx_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      end select
   end subroutine Mixed_Transform_X_Physical

   subroutine Mixed_Transform_Y_Layout(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      integer, intent(in) :: sign
      integer :: j, k, j_stop
      logical :: transform_im

      j_stop = Mixed_Active_Stop_With_Offset(1, size(phi_re_inout, 2), i_offset(npy))
      transform_im = Mixed_Input_Imag_Active(2, sign)
      select case (BC(2))
      case (ISFFT_BC_PERIODIC)
         if (2 == Mixed_First_Periodic_Axis()) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call ISFFT_1D_FFT_r2c(ny_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call ISFFT_1D_FFT_c2c(ny_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_DIRICHLET)
         if (transform_im) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call Mixed_Line_DST_Complex(ny_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call Mixed_Line_DST_Real(ny_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_NEUMANN)
         if (transform_im) then
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call Mixed_Line_DCT_Complex(ny_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, size(phi_re_inout, 3)
               do j = 1, j_stop
                  call ISFFT_1D_DCT1(ny_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      end select
   end subroutine Mixed_Transform_Y_Layout

   subroutine Mixed_Transform_Z_Layout(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      integer, intent(in) :: sign
      integer :: j, k, j_stop, k_stop
      logical :: transform_im

      j_stop = Mixed_Active_Stop_With_Offset(1, size(phi_re_inout, 2), i_offset(npy))
      k_stop = Mixed_Active_Stop_With_Offset(2, size(phi_re_inout, 3), j_offset(npz))
      transform_im = Mixed_Input_Imag_Active(3, sign)
      select case (BC(3))
      case (ISFFT_BC_PERIODIC)
         if (3 == Mixed_First_Periodic_Axis()) then
            do k = 1, k_stop
               do j = 1, j_stop
                  call ISFFT_1D_FFT_r2c(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, k_stop
               do j = 1, j_stop
                  call ISFFT_1D_FFT_c2c(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_DIRICHLET)
         if (transform_im) then
            do k = 1, k_stop
               do j = 1, j_stop
                  call Mixed_Line_DST_Complex(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, k_stop
               do j = 1, j_stop
                  call Mixed_Line_DST_Real(nz_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      case (ISFFT_BC_NEUMANN)
         if (transform_im) then
            do k = 1, k_stop
               do j = 1, j_stop
                  call Mixed_Line_DCT_Complex(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
               end do
            end do
         else
            do k = 1, k_stop
               do j = 1, j_stop
                  call ISFFT_1D_DCT1(nz_global, phi_re_inout(:,j,k), sign)
               end do
            end do
         end if
      end select
   end subroutine Mixed_Transform_Z_Layout

   subroutine Mixed_Line_DST_Real(n_global, line_re, sign)
      integer, intent(in) :: n_global, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(1:n_global)

      if (n_global > 2) call ISFFT_1D_DST1(n_global - 2, line_re(2:n_global-1), sign)
      line_re(1) = 0.0_ISFFT_REAL_KIND
      line_re(n_global) = 0.0_ISFFT_REAL_KIND
   end subroutine Mixed_Line_DST_Real

   subroutine Mixed_Line_DST_Complex(n_global, line_re, line_im, sign)
      integer, intent(in) :: n_global, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(1:n_global), line_im(1:n_global)

      if (n_global > 2) then
         call ISFFT_1D_DST1(n_global - 2, line_re(2:n_global-1), sign)
         call ISFFT_1D_DST1(n_global - 2, line_im(2:n_global-1), sign)
      end if
      line_re(1) = 0.0_ISFFT_REAL_KIND
      line_re(n_global) = 0.0_ISFFT_REAL_KIND
      line_im(1) = 0.0_ISFFT_REAL_KIND
      line_im(n_global) = 0.0_ISFFT_REAL_KIND
   end subroutine Mixed_Line_DST_Complex

   subroutine Mixed_Line_DCT_Complex(n_global, line_re, line_im, sign)
      integer, intent(in) :: n_global, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(1:n_global), line_im(1:n_global)

      call ISFFT_1D_DCT1(n_global, line_re, sign)
      call ISFFT_1D_DCT1(n_global, line_im, sign)
   end subroutine Mixed_Line_DCT_Complex

   integer function Mixed_Active_Stop_With_Offset(physical_dim, local_n, local_offset) result(stop_idx)
      integer, intent(in) :: physical_dim, local_n, local_offset
      integer :: first_periodic_axis, half_limit

      first_periodic_axis = Mixed_First_Periodic_Axis()
      if (first_periodic_axis /= physical_dim) then
         stop_idx = local_n
         return
      end if

      select case (physical_dim)
      case (1)
         half_limit = halfspectrum_x
      case (2)
         half_limit = halfspectrum_y
      case (3)
         half_limit = halfspectrum_z
      end select

      stop_idx = max(0, min(local_n, half_limit - local_offset + 1))
   end function Mixed_Active_Stop_With_Offset

   subroutine Mixed_Transform_Axis_Block(axis, bc, sign, phi_re_inout, phi_im_inout)
      integer, intent(in) :: axis, bc, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)
      integer :: ii, jj, kk, ierr
      integer :: i_stop, j_stop, k_stop
      integer :: transform_kind
      logical :: input_im_active, output_im_active, transform_im

      input_im_active = Mixed_Input_Imag_Active(axis, sign)
      output_im_active = Mixed_Output_Imag_Active(axis, sign)
      transform_im = input_im_active
      transform_kind = Mixed_Transform_Kind(bc, axis)

      select case (axis)
      case (1)
         j_stop = Mixed_Active_Local_Stop(2, axis)
         k_stop = Mixed_Active_Local_Stop(3, axis)
         do kk = 1, k_stop
            do jj = 1, j_stop
               x_sendbuf = 0.0_ISFFT_REAL_KIND
               do ii = 1, nx
                  x_sendbuf(ii + npx*nx) = phi_re_inout(ii,jj,kk)
               end do
               call MPI_Allreduce(x_sendbuf, x_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
               if (input_im_active) then
                  x_im_sendbuf = 0.0_ISFFT_REAL_KIND
                  do ii = 1, nx
                     x_im_sendbuf(ii + npx*nx) = phi_im_inout(ii,jj,kk)
                  end do
                  call MPI_Allreduce(x_im_sendbuf, x_im_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
               else
                  x_im_recvbuf = 0.0_ISFFT_REAL_KIND
               end if
               call Mixed_Transform_Dispatched_Line(transform_kind, sign, nx_global, x_recvbuf, x_im_recvbuf, transform_im)
               do ii = 1, nx
                  phi_re_inout(ii,jj,kk) = x_recvbuf(ii + npx*nx)
                  if (output_im_active) phi_im_inout(ii,jj,kk) = x_im_recvbuf(ii + npx*nx)
               end do
            end do
         end do
      case (2)
         i_stop = Mixed_Active_Local_Stop(1, axis)
         k_stop = Mixed_Active_Local_Stop(3, axis)
         do kk = 1, k_stop
            do ii = 1, i_stop
               y_sendbuf = 0.0_ISFFT_REAL_KIND
               do jj = 1, ny
                  y_sendbuf(jj + npy*ny) = phi_re_inout(ii,jj,kk)
               end do
               call MPI_Allreduce(y_sendbuf, y_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
               if (input_im_active) then
                  y_im_sendbuf = 0.0_ISFFT_REAL_KIND
                  do jj = 1, ny
                     y_im_sendbuf(jj + npy*ny) = phi_im_inout(ii,jj,kk)
                  end do
                  call MPI_Allreduce(y_im_sendbuf, y_im_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
               else
                  y_im_recvbuf = 0.0_ISFFT_REAL_KIND
               end if
               call Mixed_Transform_Dispatched_Line(transform_kind, sign, ny_global, y_recvbuf, y_im_recvbuf, transform_im)
               do jj = 1, ny
                  phi_re_inout(ii,jj,kk) = y_recvbuf(jj + npy*ny)
                  if (output_im_active) phi_im_inout(ii,jj,kk) = y_im_recvbuf(jj + npy*ny)
               end do
            end do
         end do
      case (3)
         i_stop = Mixed_Active_Local_Stop(1, axis)
         j_stop = Mixed_Active_Local_Stop(2, axis)
         do jj = 1, j_stop
            do ii = 1, i_stop
               z_sendbuf = 0.0_ISFFT_REAL_KIND
               do kk = 1, nz
                  z_sendbuf(kk + npz*nz) = phi_re_inout(ii,jj,kk)
               end do
               call MPI_Allreduce(z_sendbuf, z_recvbuf, nz_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Z, ierr)
               if (input_im_active) then
                  z_im_sendbuf = 0.0_ISFFT_REAL_KIND
                  do kk = 1, nz
                     z_im_sendbuf(kk + npz*nz) = phi_im_inout(ii,jj,kk)
                  end do
                  call MPI_Allreduce(z_im_sendbuf, z_im_recvbuf, nz_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Z, ierr)
               else
                  z_im_recvbuf = 0.0_ISFFT_REAL_KIND
               end if
               call Mixed_Transform_Dispatched_Line(transform_kind, sign, nz_global, z_recvbuf, z_im_recvbuf, transform_im)
               do kk = 1, nz
                  phi_re_inout(ii,jj,kk) = z_recvbuf(kk + npz*nz)
                  if (output_im_active) phi_im_inout(ii,jj,kk) = z_im_recvbuf(kk + npz*nz)
               end do
            end do
         end do
      end select
   end subroutine Mixed_Transform_Axis_Block

   integer function Mixed_Transform_Kind(bc, axis) result(transform_kind)
      integer, intent(in) :: bc, axis

      select case (bc)
      case (ISFFT_BC_PERIODIC)
         if (axis == Mixed_First_Periodic_Axis()) then
            transform_kind = MIXED_TRANSFORM_FFT_R2C
         else
            transform_kind = MIXED_TRANSFORM_FFT_C2C
         end if
      case (ISFFT_BC_DIRICHLET)
         transform_kind = MIXED_TRANSFORM_DST
      case (ISFFT_BC_NEUMANN)
         transform_kind = MIXED_TRANSFORM_DCT
      end select
   end function Mixed_Transform_Kind

   subroutine Mixed_Transform_Dispatched_Line(transform_kind, sign, n_global, line_re, line_im, transform_im)
      integer, intent(in) :: transform_kind, sign, n_global
      logical, intent(in) :: transform_im
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(1:n_global), line_im(1:n_global)

      select case (transform_kind)
      case (MIXED_TRANSFORM_FFT_R2C)
         call ISFFT_1D_FFT_r2c(n_global, line_re, line_im, sign)
      case (MIXED_TRANSFORM_FFT_C2C)
         call ISFFT_1D_FFT_c2c(n_global, line_re, line_im, sign)
      case (MIXED_TRANSFORM_DST)
         if (transform_im) then
            call Mixed_Line_DST_Complex(n_global, line_re, line_im, sign)
         else
            call Mixed_Line_DST_Real(n_global, line_re, sign)
         end if
      case (MIXED_TRANSFORM_DCT)
         if (transform_im) then
            call Mixed_Line_DCT_Complex(n_global, line_re, line_im, sign)
         else
            call ISFFT_1D_DCT1(n_global, line_re, sign)
         end if
      end select
   end subroutine Mixed_Transform_Dispatched_Line

   subroutine Mixed_Forward_Serial(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)
      call Mixed_Transform_Axis_Serial(1, BC(1), +1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Serial(2, BC(2), +1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Serial(3, BC(3), +1, phi_re_inout, phi_im_inout)
   end subroutine Mixed_Forward_Serial

   subroutine Mixed_Inverse_Serial(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)
      call Mixed_Transform_Axis_Serial(3, BC(3), -1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Serial(2, BC(2), -1, phi_re_inout, phi_im_inout)
      call Mixed_Transform_Axis_Serial(1, BC(1), -1, phi_re_inout, phi_im_inout)
   end subroutine Mixed_Inverse_Serial

   subroutine Mixed_Transform_Axis_Serial(axis, bc, sign, phi_re_inout, phi_im_inout)
      integer, intent(in) :: axis, bc, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(1:nx,1:ny,1:nz)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_im_inout(1:nx,1:ny,1:nz)
      integer :: i, j, k
      integer :: i_stop, j_stop, k_stop
      logical :: transform_im

      transform_im = Mixed_Input_Imag_Active(axis, sign)
      select case (axis)
      case (1)
         j_stop = Mixed_Active_Local_Stop(2, axis)
         k_stop = Mixed_Active_Local_Stop(3, axis)
         select case (bc)
         case (ISFFT_BC_PERIODIC)
            if (axis == Mixed_First_Periodic_Axis()) then
               do k = 1, k_stop
                  do j = 1, j_stop
                     call ISFFT_1D_FFT_r2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do j = 1, j_stop
                     call ISFFT_1D_FFT_c2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_DIRICHLET)
            if (transform_im) then
               do k = 1, k_stop
                  do j = 1, j_stop
                     call Mixed_Line_DST_Complex(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do j = 1, j_stop
                     call Mixed_Line_DST_Real(nx_global, phi_re_inout(:,j,k), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_NEUMANN)
            if (transform_im) then
               do k = 1, k_stop
                  do j = 1, j_stop
                     call Mixed_Line_DCT_Complex(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do j = 1, j_stop
                     call ISFFT_1D_DCT1(nx_global, phi_re_inout(:,j,k), sign)
                  end do
               end do
            end if
         end select
      case (2)
         i_stop = Mixed_Active_Local_Stop(1, axis)
         k_stop = Mixed_Active_Local_Stop(3, axis)
         select case (bc)
         case (ISFFT_BC_PERIODIC)
            if (axis == Mixed_First_Periodic_Axis()) then
               do k = 1, k_stop
                  do i = 1, i_stop
                     call ISFFT_1D_FFT_r2c(ny_global, phi_re_inout(i,:,k), phi_im_inout(i,:,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do i = 1, i_stop
                     call ISFFT_1D_FFT_c2c(ny_global, phi_re_inout(i,:,k), phi_im_inout(i,:,k), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_DIRICHLET)
            if (transform_im) then
               do k = 1, k_stop
                  do i = 1, i_stop
                     call Mixed_Line_DST_Complex(ny_global, phi_re_inout(i,:,k), phi_im_inout(i,:,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do i = 1, i_stop
                     call Mixed_Line_DST_Real(ny_global, phi_re_inout(i,:,k), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_NEUMANN)
            if (transform_im) then
               do k = 1, k_stop
                  do i = 1, i_stop
                     call Mixed_Line_DCT_Complex(ny_global, phi_re_inout(i,:,k), phi_im_inout(i,:,k), sign)
                  end do
               end do
            else
               do k = 1, k_stop
                  do i = 1, i_stop
                     call ISFFT_1D_DCT1(ny_global, phi_re_inout(i,:,k), sign)
                  end do
               end do
            end if
         end select
      case (3)
         i_stop = Mixed_Active_Local_Stop(1, axis)
         j_stop = Mixed_Active_Local_Stop(2, axis)
         select case (bc)
         case (ISFFT_BC_PERIODIC)
            if (axis == Mixed_First_Periodic_Axis()) then
               do j = 1, j_stop
                  do i = 1, i_stop
                     call ISFFT_1D_FFT_r2c(nz_global, phi_re_inout(i,j,:), phi_im_inout(i,j,:), sign)
                  end do
               end do
            else
               do j = 1, j_stop
                  do i = 1, i_stop
                     call ISFFT_1D_FFT_c2c(nz_global, phi_re_inout(i,j,:), phi_im_inout(i,j,:), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_DIRICHLET)
            if (transform_im) then
               do j = 1, j_stop
                  do i = 1, i_stop
                     call Mixed_Line_DST_Complex(nz_global, phi_re_inout(i,j,:), phi_im_inout(i,j,:), sign)
                  end do
               end do
            else
               do j = 1, j_stop
                  do i = 1, i_stop
                     call Mixed_Line_DST_Real(nz_global, phi_re_inout(i,j,:), sign)
                  end do
               end do
            end if
         case (ISFFT_BC_NEUMANN)
            if (transform_im) then
               do j = 1, j_stop
                  do i = 1, i_stop
                     call Mixed_Line_DCT_Complex(nz_global, phi_re_inout(i,j,:), phi_im_inout(i,j,:), sign)
                  end do
               end do
            else
               do j = 1, j_stop
                  do i = 1, i_stop
                     call ISFFT_1D_DCT1(nz_global, phi_re_inout(i,j,:), sign)
                  end do
               end do
            end if
         end select
      end select
   end subroutine Mixed_Transform_Axis_Serial

   integer function Mixed_First_Periodic_Axis() result(axis)
      integer :: d

      axis = 0
      do d = 1, 3
         if (BC(d) == ISFFT_BC_PERIODIC) then
            axis = d
            return
         end if
      end do
   end function Mixed_First_Periodic_Axis

   integer function Mixed_Active_Local_Stop(dim, transform_axis) result(stop_idx)
      integer, intent(in) :: dim, transform_axis
      integer :: first_periodic_axis, half_limit, local_offset, local_n

      select case (dim)
      case (1)
         local_n = nx
         local_offset = i_offset(npx)
         half_limit = halfspectrum_x
      case (2)
         local_n = ny
         local_offset = j_offset(npy)
         half_limit = halfspectrum_y
      case (3)
         local_n = nz
         local_offset = k_offset(npz)
         half_limit = halfspectrum_z
      end select

      first_periodic_axis = Mixed_First_Periodic_Axis()
      if (first_periodic_axis == dim .and. transform_axis > first_periodic_axis) then
         stop_idx = max(0, min(local_n, half_limit - local_offset + 1))
      else
         stop_idx = local_n
      end if
   end function Mixed_Active_Local_Stop

   logical function Mixed_Any_Periodic() result(has_periodic)
      has_periodic = any(BC == ISFFT_BC_PERIODIC)
   end function Mixed_Any_Periodic

   logical function Mixed_Input_Imag_Active(axis, sign) result(active)
      integer, intent(in) :: axis, sign

      if (sign > 0) then
         active = any(BC(1:axis-1) == ISFFT_BC_PERIODIC)
      else
         active = any(BC(1:axis) == ISFFT_BC_PERIODIC)
      end if
   end function Mixed_Input_Imag_Active

   logical function Mixed_Output_Imag_Active(axis, sign) result(active)
      integer, intent(in) :: axis, sign

      if (sign > 0) then
         active = any(BC(1:axis) == ISFFT_BC_PERIODIC)
      else if (axis > 1) then
         active = any(BC(1:axis-1) == ISFFT_BC_PERIODIC)
      else
         active = .false.
      end if
   end function Mixed_Output_Imag_Active

   subroutine Mixed_Algebra_Block(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      type(ISFFT_Scheme_Coefficients) :: coeff
      integer :: kx, ky, kz, k1, k2, k3
      logical :: on_dirichlet_boundary, on_x_boundary, on_y_boundary, on_z_boundary
      real(kind=ISFFT_REAL_KIND) :: lambda_x, lambda_y, lambda_z, lambda_total
      real(kind=ISFFT_REAL_KIND), allocatable :: lambda_xs(:), lambda_ys(:), lambda_zs(:)
      logical, allocatable :: bx(:), by(:), bz(:), ax(:), ay(:), az(:)

      call ISFFT_Get_Scheme_Coefficients(If_scheme, coeff)

      kx = size(phi_re_inout, 1)
      ky = size(phi_re_inout, 2)
      kz = size(phi_re_inout, 3)

      allocate(lambda_xs(kx), lambda_ys(ky), lambda_zs(kz))
      allocate(bx(kx), by(ky), bz(kz), ax(kx), ay(ky), az(kz))
      call Mixed_Precompute_Mode_Info(1, kx, npx*kx, coeff, lambda_xs, bx, ax)
      call Mixed_Precompute_Mode_Info(2, ky, npy*ky, coeff, lambda_ys, by, ay)
      call Mixed_Precompute_Mode_Info(3, kz, npz*kz, coeff, lambda_zs, bz, az)

      do k3 = 1, kz
         if (.not. az(k3)) cycle
         on_z_boundary = bz(k3)
         lambda_z = lambda_zs(k3)
         do k2 = 1, ky
            if (.not. ay(k2)) cycle
            on_y_boundary = by(k2)
            lambda_y = lambda_ys(k2)
            do k1 = 1, kx
               if (.not. ax(k1)) cycle
               on_x_boundary = bx(k1)
               lambda_x = lambda_xs(k1)

               on_dirichlet_boundary = on_x_boundary .or. on_y_boundary .or. on_z_boundary
               lambda_total = lambda_x + lambda_y + lambda_z

               if (.not. on_dirichlet_boundary .and. lambda_total /= 0.0_ISFFT_REAL_KIND) then
                  phi_re_inout(k1,k2,k3) = -phi_re_inout(k1,k2,k3)/lambda_total
                  if (Mixed_Any_Periodic()) phi_im_inout(k1,k2,k3) = -phi_im_inout(k1,k2,k3)/lambda_total
               else
                  phi_re_inout(k1,k2,k3) = 0.0_ISFFT_REAL_KIND
                  if (Mixed_Any_Periodic()) phi_im_inout(k1,k2,k3) = 0.0_ISFFT_REAL_KIND
               end if
            end do
         end do
      end do

      deallocate(lambda_xs, lambda_ys, lambda_zs, bx, by, bz, ax, ay, az)
   end subroutine Mixed_Algebra_Block

   subroutine Mixed_Algebra_ZXY(phi_re_inout, phi_im_inout)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      type(ISFFT_Scheme_Coefficients) :: coeff
      integer :: kz, kx, ky, k1, k2, k3
      logical :: on_dirichlet_boundary, on_x_boundary, on_y_boundary, on_z_boundary
      real(kind=ISFFT_REAL_KIND) :: lambda_x, lambda_y, lambda_z, lambda_total
      real(kind=ISFFT_REAL_KIND), allocatable :: lambda_xs(:), lambda_ys(:), lambda_zs(:)
      logical, allocatable :: bx(:), by(:), bz(:), ax(:), ay(:), az(:)

      call ISFFT_Get_Scheme_Coefficients(If_scheme, coeff)

      kz = size(phi_re_inout, 1)
      kx = size(phi_re_inout, 2)
      ky = size(phi_re_inout, 3)

      allocate(lambda_zs(kz), lambda_xs(kx), lambda_ys(ky))
      allocate(bz(kz), bx(kx), by(ky), az(kz), ax(kx), ay(ky))
      call Mixed_Precompute_Mode_Info(3, kz, npx*kz, coeff, lambda_zs, bz, az)
      call Mixed_Precompute_Mode_Info(1, kx, npy*kx, coeff, lambda_xs, bx, ax)
      call Mixed_Precompute_Mode_Info(2, ky, npz*ky, coeff, lambda_ys, by, ay)

      do k3 = 1, ky
         if (.not. ay(k3)) cycle
         on_y_boundary = by(k3)
         lambda_y = lambda_ys(k3)
         do k2 = 1, kx
            if (.not. ax(k2)) cycle
            on_x_boundary = bx(k2)
            lambda_x = lambda_xs(k2)
            do k1 = 1, kz
               if (.not. az(k1)) cycle
               on_z_boundary = bz(k1)
               lambda_z = lambda_zs(k1)

               on_dirichlet_boundary = on_x_boundary .or. on_y_boundary .or. on_z_boundary
               lambda_total = lambda_x + lambda_y + lambda_z

               if (.not. on_dirichlet_boundary .and. lambda_total /= 0.0_ISFFT_REAL_KIND) then
                  phi_re_inout(k1,k2,k3) = -phi_re_inout(k1,k2,k3)/lambda_total
                  if (Mixed_Any_Periodic()) phi_im_inout(k1,k2,k3) = -phi_im_inout(k1,k2,k3)/lambda_total
               else
                  phi_re_inout(k1,k2,k3) = 0.0_ISFFT_REAL_KIND
                  if (Mixed_Any_Periodic()) phi_im_inout(k1,k2,k3) = 0.0_ISFFT_REAL_KIND
               end if
            end do
         end do
      end do

      deallocate(lambda_zs, lambda_xs, lambda_ys, bz, bx, by, az, ax, ay)
   end subroutine Mixed_Algebra_ZXY

   subroutine Mixed_Precompute_Mode_Info(dim, local_n, local_offset, coeff, lambda, on_boundary, active)
      integer, intent(in) :: dim, local_n, local_offset
      type(ISFFT_Scheme_Coefficients), intent(in) :: coeff
      real(kind=ISFFT_REAL_KIND), intent(out) :: lambda(local_n)
      logical, intent(out) :: on_boundary(local_n), active(local_n)
      integer :: idx, global_idx, mode_idx, n_global
      real(kind=ISFFT_REAL_KIND) :: length

      select case (dim)
      case (1)
         n_global = nx_global
         length = SLx
      case (2)
         n_global = ny_global
         length = SLy
      case (3)
         n_global = nz_global
         length = SLz
      end select

      do idx = 1, local_n
         global_idx = idx + local_offset
         active(idx) = Mixed_Global_Mode_Active(dim, global_idx)
         call Mixed_Mode_Info(BC(dim), global_idx, n_global, mode_idx, on_boundary(idx))
         lambda(idx) = ISFFT_Lambda_For_BC(BC(dim), mode_idx, n_global, length, If_scheme, coeff)
      end do
   end subroutine Mixed_Precompute_Mode_Info

   logical function Mixed_Global_Mode_Active(dim, global_index_one_based) result(active)
      integer, intent(in) :: dim, global_index_one_based
      integer :: first_periodic_axis, half_limit

      first_periodic_axis = Mixed_First_Periodic_Axis()
      if (first_periodic_axis /= dim) then
         active = .true.
         return
      end if

      select case (dim)
      case (1)
         half_limit = halfspectrum_x
      case (2)
         half_limit = halfspectrum_y
      case (3)
         half_limit = halfspectrum_z
      end select
      active = (global_index_one_based <= half_limit)
   end function Mixed_Global_Mode_Active

   subroutine Mixed_Mode_Info(bc, global_index_one_based, n_global, mode_index, on_dirichlet_boundary)
      integer, intent(in) :: bc, global_index_one_based, n_global
      integer, intent(out) :: mode_index
      logical, intent(out) :: on_dirichlet_boundary
      integer :: raw_mode

      raw_mode = global_index_one_based - 1
      on_dirichlet_boundary = .false.

      select case (bc)
      case (ISFFT_BC_PERIODIC)
         mode_index = raw_mode
         if (2*mode_index >= n_global) mode_index = mode_index - n_global
      case (ISFFT_BC_DIRICHLET)
         on_dirichlet_boundary = (global_index_one_based == 1 .or. global_index_one_based == n_global)
         mode_index = merge(0, raw_mode, on_dirichlet_boundary)
      case (ISFFT_BC_NEUMANN)
         mode_index = raw_mode
      case default
         mode_index = raw_mode
      end select
   end subroutine Mixed_Mode_Info

   subroutine Mixed_Dirichlet_HomoBC(phi)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      integer :: i, j, k, gi, gj, gk

      if (BC(1) == ISFFT_BC_DIRICHLET) then
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  gi = i_offset(npx) + i - 1
                  if (gi == 1 .or. gi == nx_global) phi(i,j,k) = 0.0_ISFFT_REAL_KIND
               end do
            end do
         end do
      end if

      if (BC(2) == ISFFT_BC_DIRICHLET) then
         do k = 1, nz
            do j = 1, ny
               gj = j_offset(npy) + j - 1
               if (gj == 1 .or. gj == ny_global) phi(:,j,k) = 0.0_ISFFT_REAL_KIND
            end do
         end do
      end if

      if (BC(3) == ISFFT_BC_DIRICHLET) then
         do k = 1, nz
            gk = k_offset(npz) + k - 1
            if (gk == 1 .or. gk == nz_global) phi(:,:,k) = 0.0_ISFFT_REAL_KIND
         end do
      end if
   end subroutine Mixed_Dirichlet_HomoBC

end module ISFFT_3D_Mixed_Poisson
