module SPLASH_3D_FFT
   use SPLASH_Parameters
   use SPLASH_MPI_Constants
   use SPLASH_Buffer
   use SPLASH_LocalTranspose
   use SPLASH_GlobalTranspose
   use SPLASH_FFT_kernel
   implicit none

contains

   subroutine FFT_x_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------xFFT + MPI_allreduce comm_x------------------------------
      x_sendbuf = 0.d0
      x_im_sendbuf = 0.d0
      x_recvbuf = 0.d0
      x_im_recvbuf = 0.d0

      do kk = 1, nz
         do jj = 1, ny

            do ii = 1, nx
               x_sendbuf(ii + npx*nx) = phi_re_inout(ii, jj, kk)
               x_im_sendbuf(ii + npx*nx) = phi_im_inout(ii, jj, kk)
            enddo

            call MPI_allREDUCE(x_sendbuf, x_recvbuf, nx_global, SPLASH_DATA_TYPE, MPI_SUM, &
               MPI_COMM_X, ierr)
            call MPI_allREDUCE(x_im_sendbuf, x_im_recvbuf, nx_global, SPLASH_DATA_TYPE, MPI_SUM, &
               MPI_COMM_X, ierr)

            call SPLASH_1D_FFT_r2c(nx_global, x_recvbuf, x_im_recvbuf, 1)

            do ii = 1, nx
               phi_re_inout(ii, jj, kk) = x_recvbuf(ii + npx*nx)
               phi_im_inout(ii, jj, kk) = x_im_recvbuf(ii + npx*nx)
            enddo

         enddo
      enddo

   end subroutine FFT_x_block
   !----------------------------------------------------------------------
   subroutine FFT_y_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------yFFT + MPI_allreduce comm_y------------------------------!
      y_sendbuf = 0.d0
      y_im_sendbuf = 0.d0
      y_recvbuf = 0.d0
      y_im_recvbuf = 0.d0

      if (i_offset(npx) .le. halfspectrum_x) then
         do kk = 1, nz
            do ii = 1, min(nx,(halfspectrum_x - i_offset(npx) + 1))

               do jj = 1, ny
                  y_sendbuf(jj + npy*ny) = phi_re_inout(ii, jj, kk)
                  y_im_sendbuf(jj + npy*ny) = phi_im_inout(ii, jj, kk)
               enddo

               call MPI_allREDUCE(y_sendbuf, y_recvbuf, ny_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Y, ierr)
               call MPI_allREDUCE(y_im_sendbuf, y_im_recvbuf, ny_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Y, ierr)

               call SPLASH_1D_FFT_c2c(ny_global, y_recvbuf, y_im_recvbuf, 1)

               do jj = 1, ny
                  phi_re_inout(ii, jj, kk) = y_recvbuf(jj + npy*ny)
                  phi_im_inout(ii, jj, kk) = y_im_recvbuf(jj + npy*ny)
               enddo

            enddo
         enddo
      endif

   end subroutine FFT_y_block
   !----------------------------------------------------------------------
   subroutine FFT_z_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------zFFT + MPI_allreduce comm_z------------------------------!
      z_sendbuf = 0.d0
      z_im_sendbuf = 0.d0
      z_recvbuf = 0.d0
      z_im_recvbuf = 0.d0

      if (i_offset(npx) .le. halfspectrum_x) then
         do jj = 1, ny
            do ii = 1, min(nx,(halfspectrum_x - i_offset(npx) + 1))

               do kk = 1, nz
                  z_sendbuf(kk + npz*nz) = phi_re_inout(ii, jj, kk)
                  z_im_sendbuf(kk + npz*nz) = phi_im_inout(ii, jj, kk)
               enddo

               call MPI_allREDUCE(z_sendbuf, z_recvbuf, nz_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Z, ierr)
               call MPI_allREDUCE(z_im_sendbuf, z_im_recvbuf, nz_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Z, ierr)

               call SPLASH_1D_FFT_c2c(nz_global, z_recvbuf, z_im_recvbuf, 1)

               do kk = 1, nz
                  phi_re_inout(ii, jj, kk) = z_recvbuf(kk + npz*nz)
                  phi_im_inout(ii, jj, kk) = z_im_recvbuf(kk + npz*nz)
               enddo

            enddo
         enddo
      endif

   end subroutine FFT_z_block
   !----------------------------------------------------------------------
   subroutine iFFT_z_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------zIFFT + MPI_allreduce comm_z------------------------------!
      z_sendbuf = 0.d0
      z_im_sendbuf = 0.d0
      z_recvbuf = 0.d0
      z_im_recvbuf = 0.d0

      if (i_offset(npx) .le. halfspectrum_x) then
         do jj = 1, ny
            do ii = 1, min(nx,(halfspectrum_x - i_offset(npx) + 1))

               do kk = 1, nz
                  z_sendbuf(kk + npz*nz) = phi_re_inout(ii, jj, kk)
                  z_im_sendbuf(kk + npz*nz) = phi_im_inout(ii, jj, kk)
               enddo

               call MPI_allREDUCE(z_sendbuf, z_recvbuf, nz_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Z, ierr)
               call MPI_allREDUCE(z_im_sendbuf, z_im_recvbuf, nz_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Z, ierr)

               call SPLASH_1D_FFT_c2c(nz_global, z_recvbuf, z_im_recvbuf, -1)

               do kk = 1, nz
                  phi_re_inout(ii, jj, kk) = z_recvbuf(kk + npz*nz)
                  phi_im_inout(ii, jj, kk) = z_im_recvbuf(kk + npz*nz)
               enddo

            enddo
         enddo
      endif

   end subroutine iFFT_z_block
   !----------------------------------------------------------------------
   subroutine iFFT_y_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------yIFFT + MPI_allreduce comm_y------------------------------!
      y_sendbuf = 0.d0
      y_im_sendbuf = 0.d0
      y_recvbuf = 0.d0
      y_im_recvbuf = 0.d0

      if (i_offset(npx) .le. halfspectrum_x) then
         do kk = 1, nz
            do ii = 1, min(nx,(halfspectrum_x - i_offset(npx) + 1))

               do jj = 1, ny
                  y_sendbuf(jj + npy*ny) = phi_re_inout(ii, jj, kk)
                  y_im_sendbuf(jj + npy*ny) = phi_im_inout(ii, jj, kk)
               enddo

               call MPI_allREDUCE(y_sendbuf, y_recvbuf, ny_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Y, ierr)
               call MPI_allREDUCE(y_im_sendbuf, y_im_recvbuf, ny_global, SPLASH_DATA_TYPE, MPI_SUM, &
                  MPI_COMM_Y, ierr)

               call SPLASH_1D_FFT_c2c(ny_global, y_recvbuf, y_im_recvbuf, -1)

               do jj = 1, ny
                  phi_re_inout(ii, jj, kk) = y_recvbuf(jj + npy*ny)
                  phi_im_inout(ii, jj, kk) = y_im_recvbuf(jj + npy*ny)
               enddo

            enddo
         enddo
      endif

   end subroutine iFFT_y_block
   !----------------------------------------------------------------------
   subroutine iFFT_x_block(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      integer :: ii, jj, kk, ierr
      !-------------------------xIFFT + MPI_allreduce comm_x------------------------------!
      x_sendbuf = 0.d0
      x_im_sendbuf = 0.d0
      x_recvbuf = 0.d0
      x_im_recvbuf = 0.d0

      do kk = 1, nz
         do jj = 1, ny

            do ii = 1, nx
               x_sendbuf(ii + npx*nx) = phi_re_inout(ii, jj, kk)
               x_im_sendbuf(ii + npx*nx) = phi_im_inout(ii, jj, kk)
            enddo

            call MPI_allREDUCE(x_sendbuf, x_recvbuf, nx_global, SPLASH_DATA_TYPE, MPI_SUM, &
               MPI_COMM_X, ierr)
            call MPI_allREDUCE(x_im_sendbuf, x_im_recvbuf, nx_global, SPLASH_DATA_TYPE, MPI_SUM, &
               MPI_COMM_X, ierr)

            call SPLASH_1D_FFT_r2c(nx_global, x_recvbuf, x_im_recvbuf, -1)

            do ii = 1, nx
               phi_re_inout(ii, jj, kk) = x_recvbuf(ii + npx*nx)
               phi_im_inout(ii, jj, kk) = x_im_recvbuf(ii + npx*nx)
            enddo

         enddo
      enddo

   end subroutine iFFT_x_block
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine FFT_x_slab(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      ! real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------xFFT + 0*lt + 0*gt------------------------------!
      do k = 1, nz
         do j = 1, ny
            call SPLASH_1D_FFT_r2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), 1)
         enddo
      enddo

   end subroutine FFT_x_slab
   !----------------------------------------------------------------------
   subroutine FFT_y_slab(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------yFFT + 1*lt + 0*gt------------------------------!
      call local_transpose_12(phi_re_in,phi_im_in,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, nz
            do j = 1, min(nx,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(ny_global, phi_re_out(:,j,k), phi_im_out(:,j,k), 1)
            enddo
         enddo
      endif

   end subroutine FFT_y_slab
   !----------------------------------------------------------------------
   subroutine FFT_z_slab(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:ny, 1:nx, 1:nz) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------zFFT + 1*lt + 1*gt------------------------------!
      call local_transpose_13(phi_re_in,phi_im_in,phi_re_transposed_a,phi_im_transposed_a)

      call global_transpose_13(phi_re_transposed_a,phi_im_transposed_a,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, ny_global/npz0
            do j = 1, min(nx_global/npx0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(nz_global, phi_re_out(:,j,k), phi_im_out(:,j,k), 1)
            enddo
         enddo
      endif

   end subroutine FFT_z_slab
   !----------------------------------------------------------------------
   subroutine iFFT_z_slab(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nz_global/npy0, 1:nx_global/npx0, 1:ny_global/npz0) &
         , intent(inout) :: phi_re_inout, phi_im_inout
      ! real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------zIFFT + 0*lt + 0*gt------------------------------!
      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, ny_global/npz0
            do j = 1, min(nx_global/npx0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), -1)
            enddo
         enddo
      endif

   end subroutine iFFT_z_slab
   !----------------------------------------------------------------------
   subroutine iFFT_y_slab(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:nz_global/npy0, 1:nx_global/npx0, 1:ny_global/npz0) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------yIFFT + 1*lt + 1*gt------------------------------!
      call global_transpose_31(phi_re_in,phi_im_in,phi_re_transposed_a,phi_im_transposed_a)

      call local_transpose_13(phi_re_transposed_a,phi_im_transposed_a,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, nz_global/npz0
            do j = 1, min(nx_global/npx0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(ny_global, phi_re_out(:,j,k), phi_im_out(:,j,k), -1)
            enddo
         enddo
      endif

   end subroutine iFFT_y_slab
   !----------------------------------------------------------------------
   subroutine iFFT_x_slab(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:ny_global/npy0, 1:nx_global/npx0, 1:nz_global/npz0) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out

      integer :: i, j, k, ierr
      !-------------------------xIFFT + 1*lt + 0*gt------------------------------!
      call local_transpose_12(phi_re_in,phi_im_in,phi_re_out,phi_im_out)

      do k = 1, nz_global/npz0
         do j = 1, ny_global/npy0
            call SPLASH_1D_FFT_r2c(nx_global, phi_re_out(:,j,k), phi_im_out(:,j,k), -1)
         enddo
      enddo

   end subroutine iFFT_x_slab
   !----------------------------------------------------------------------


   !----------------------------------------------------------------------
   subroutine FFT_x_pencil(phi_re_inout, phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz), intent(inout) :: phi_re_inout, phi_im_inout
      ! real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------xFFT + 0*lt + 0*gt------------------------------!
      do k = 1, nz
         do j = 1, ny
            call SPLASH_1D_FFT_r2c(nx_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), 1)
         enddo
      enddo

   end subroutine FFT_x_pencil
   !----------------------------------------------------------------------
   subroutine FFT_y_pencil(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:nx, 1:ny, 1:nz) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------yFFT + 1*lt + 1*gt------------------------------!
      call local_transpose_12(phi_re_in,phi_im_in,phi_re_transposed_a,phi_im_transposed_a)

      call global_transpose_12(phi_re_transposed_a,phi_im_transposed_a,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, nz_global/npz0
            do j = 1, min(nx_global/npy0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(ny_global, phi_re_out(:,j,k), phi_im_out(:,j,k), 1)
            enddo
         enddo
      endif

   end subroutine FFT_y_pencil
   !----------------------------------------------------------------------
   subroutine FFT_z_pencil(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:ny_global/npx0, 1:nx_global/npy0, 1:nz_global/npz0) &
         , intent(inout) :: phi_re_in,phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------zFFT + 1*lt + 1*gt------------------------------!
      call local_transpose_13(phi_re_in,phi_im_in,phi_re_transposed_b,phi_im_transposed_b)

      call global_transpose_13(phi_re_transposed_b,phi_im_transposed_b,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, ny_global/npz0
            do j = 1, min(nx_global/npy0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(nz_global, phi_re_out(:,j,k), phi_im_out(:,j,k), 1)
            enddo
         enddo
      endif

   end subroutine FFT_z_pencil
   !----------------------------------------------------------------------
   subroutine iFFT_z_pencil(phi_re_inout,phi_im_inout)
      real(kind=SPLASH_REAL_KIND), dimension(1:nz_global/npx0, 1:nx_global/npy0, 1:ny_global/npz0) &
         , intent(inout) :: phi_re_inout, phi_im_inout
      ! real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------zIFFT + 0*lt + 0*gt------------------------------!
      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, ny_global/npz0
            do j = 1, min(nx_global/npy0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(nz_global, phi_re_inout(:,j,k), phi_im_inout(:,j,k), -1)
            enddo
         enddo
      endif

   end subroutine iFFT_z_pencil
   !----------------------------------------------------------------------
   subroutine iFFT_y_pencil(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:nz_global/npx0, 1:nx_global/npy0, 1:ny_global/npz0) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------yIFFT + 1*lt + 1*gt------------------------------!
      call global_transpose_31(phi_re_in,phi_im_in,phi_re_transposed_b,phi_im_transposed_b)

      call local_transpose_13(phi_re_transposed_b,phi_im_transposed_b,phi_re_out,phi_im_out)

      if (i_offset(npx) .le. halfspectrum_x) then
         do k = 1, nz_global/npz0
            do j = 1, min(nx_global/npy0,(halfspectrum_x - i_offset(npx) + 1))
               call SPLASH_1D_FFT_c2c(ny_global, phi_re_out(:,j,k), phi_im_out(:,j,k), -1)
            enddo
         enddo
      endif

   end subroutine iFFT_y_pencil
   !----------------------------------------------------------------------
   subroutine iFFT_x_pencil(phi_re_in,phi_im_in,phi_re_out,phi_im_out)
      real(kind=SPLASH_REAL_KIND), dimension(1:ny_global/npx0, 1:nx_global/npy0, 1:nz_global/npz0) &
         , intent(inout) :: phi_re_in, phi_im_in
      real(kind=SPLASH_REAL_KIND), dimension(:,:,:), intent(inout) :: phi_re_out, phi_im_out
      integer :: i, j, k, ierr
      !-------------------------xIFFT + 1*lt + 1*gt------------------------------!
      call global_transpose_21(phi_re_in,phi_im_in,phi_re_transposed_a,phi_im_transposed_a)

      call local_transpose_12(phi_re_transposed_a,phi_im_transposed_a,phi_re_out,phi_im_out)

      do k = 1, nz_global/npz0
         do j = 1, ny_global/npy0
            call SPLASH_1D_FFT_r2c(nx_global, phi_re_out(:,j,k), phi_im_out(:,j,k), -1)
         enddo
      enddo

   end subroutine iFFT_x_pencil
   !----------------------------------------------------------------------





end module SPLASH_3D_FFT
