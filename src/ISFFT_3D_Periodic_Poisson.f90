!!----------------------------------------------------------------------
module ISFFT_3D_Periodic_Poisson
   use ISFFT_Buffer
   use ISFFT_3D_FFT
   use ISFFT_Poisson_Algebra
   implicit none

   private
   public :: ISFFT_3D_Periodic_Poisson_serial, &
      ISFFT_3D_Periodic_Poisson_block, &
      ISFFT_3D_Periodic_Poisson_slab, &
      ISFFT_3D_Periodic_Poisson_pencil

contains

!!---------------------------For test only------------------------------
   Subroutine ISFFT_3D_Periodic_Poisson_serial(phi,f)
      implicit none
      integer :: i, j, k

      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: f
      real(kind=ISFFT_REAL_KIND), dimension(1:nx,1:ny,1:nz) :: f_re_serial, f_im_serial
      !!----------------------------------------------------------------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f_re_serial(i,j,k)=f(i,j,k)
            enddo
         enddo
      enddo

      do k = 1, nz
         do j = 1, ny
            call ISFFT_1D_FFT_r2c(nx_global, f_re_serial(:,j,k), f_im_serial(:,j,k), 1)
         enddo
      enddo

      do k = 1, nz
         do j = 1, nx/2 + 1
            call ISFFT_1D_FFT_c2c(ny_global, f_re_serial(j,:,k), f_im_serial(j,:,k), 1)
         enddo
      enddo

      do k = 1, ny
         do j = 1, nx/2 + 1
            call ISFFT_1D_FFT_c2c(nz_global, f_re_serial(j,k,:), f_im_serial(j,k,:), 1)
         enddo
      enddo

      call Algebra_Periodic_Poisson_r2c_serial(f_re_serial,f_im_serial)

      do k = 1, ny
         do j = 1, nx/2 + 1
            call ISFFT_1D_FFT_c2c(nz_global, f_re_serial(j,k,:), f_im_serial(j,k,:), -1)
         enddo
      enddo

      do k = 1, nz
         do j = 1, nx/2 + 1
            call ISFFT_1D_FFT_c2c(ny_global, f_re_serial(j,:,k), f_im_serial(j,:,k), -1)
         enddo
      enddo

      do k = 1, nz
         do j = 1, ny
            call ISFFT_1D_FFT_r2c(nx_global, f_re_serial(:,j,k), f_im_serial(:,j,k), -1)
         enddo
      enddo

      do k=1,nz
         do j=1,ny
            do i=1,nx
               phi(i,j,k)=f_re_serial(i,j,k)
            enddo
         enddo
      enddo

      if (my_id .eq. 0) print*, 'Fourier_one_3D_serial is done'

   End Subroutine ISFFT_3D_Periodic_Poisson_serial
!!----------------------------------------------------------------------
   Subroutine ISFFT_3D_Periodic_Poisson_block(phi,f)
      implicit none
      integer :: i, j, k, ierr

      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: f
      !!----------------------------------------------------------------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f_re(i,j,k)=f(i,j,k)
            enddo
         enddo
      enddo

      call FFT_x_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_y_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_z_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Algebra_Periodic_Poisson_r2c(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_z_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_y_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_x_block(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               phi(i,j,k)=f_re(i,j,k)
            enddo
         enddo
      enddo

   End Subroutine ISFFT_3D_Periodic_Poisson_block
!!----------------------------------------------------------------------
   Subroutine ISFFT_3D_Periodic_Poisson_slab(phi,f)
      implicit none
      integer :: i, j, k, ierr
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: f

      !!----------------------------------------------------------------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f_re(i,j,k)=f(i,j,k)
            enddo
         enddo
      enddo

      call FFT_x_slab(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_y_slab(f_re, f_im, f_hat_y,  f_hat_im_y)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_z_slab(f_hat_y, f_hat_im_y, f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Algebra_Periodic_Poisson_r2c(f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_z_slab(f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_y_slab(f_hat, f_hat_im, phi_hat_x, phi_hat_im_x)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_x_slab(phi_hat_x, phi_hat_im_x, phi_re, phi_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               phi(i,j,k)=phi_re(i,j,k)
            enddo
         enddo
      enddo

   End Subroutine ISFFT_3D_Periodic_Poisson_slab
!!----------------------------------------------------------------------
   Subroutine ISFFT_3D_Periodic_Poisson_pencil(phi,f)
      implicit none
      integer :: i, j, k, ierr

      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: phi
      real(kind=ISFFT_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: f
      !!----------------------------------------------------------------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f_re(i,j,k)=f(i,j,k)
            enddo
         enddo
      enddo

      call FFT_x_pencil(f_re, f_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_y_pencil(f_re, f_im, f_hat_y, f_hat_im_y)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call FFT_z_pencil(f_hat_y, f_hat_im_y, f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call Algebra_Periodic_Poisson_r2c(f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_z_pencil(f_hat, f_hat_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_y_pencil(f_hat, f_hat_im, phi_hat_x, phi_hat_im_x)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call iFFT_x_pencil(phi_hat_x, phi_hat_im_x, phi_re, phi_im)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               phi(i,j,k)=phi_re(i,j,k)
            enddo
         enddo
      enddo

   End Subroutine ISFFT_3D_Periodic_Poisson_pencil
!!----------------------------------------------------------------------


end module ISFFT_3D_Periodic_Poisson
!!----------------------------------------------------------------------
