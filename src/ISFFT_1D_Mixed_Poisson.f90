module ISFFT_1D_Mixed_Poisson
   use ISFFT_Parameters
   use ISFFT_FFT_kernel, only: ISFFT_1D_FFT_r2c
   use ISFFT_DST_kernel, only: ISFFT_1D_DST1
   use ISFFT_DCT_kernel, only: ISFFT_1D_DCT1
   use ISFFT_Discrete_Symbols
   implicit none

   private
   public :: ISFFT_1D_Mixed_Poisson_serial

contains

   subroutine ISFFT_1D_Mixed_Poisson_serial(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: &
         phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(in) :: &
         f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND) :: work_re(nx_global), work_im(nx_global)
      integer :: i

      work_re = f(1:nx_global,1,1)
      work_im = 0.0_ISFFT_REAL_KIND

      select case (BC(1))
      case (ISFFT_BC_PERIODIC)
         call ISFFT_1D_FFT_r2c(nx_global, work_re, work_im, +1)
         call Algebra_Periodic(work_re, work_im)
         call ISFFT_1D_FFT_r2c(nx_global, work_re, work_im, -1)
      case (ISFFT_BC_DIRICHLET)
         work_re(1) = 0.0_ISFFT_REAL_KIND
         work_re(nx_global) = 0.0_ISFFT_REAL_KIND
         if (nx_global > 2) call ISFFT_1D_DST1(nx_global-2, work_re(2:nx_global-1), +1)
         call Algebra_Real(work_re)
         if (nx_global > 2) call ISFFT_1D_DST1(nx_global-2, work_re(2:nx_global-1), -1)
      case (ISFFT_BC_NEUMANN)
         call ISFFT_1D_DCT1(nx_global, work_re, +1)
         call Algebra_Real(work_re)
         call ISFFT_1D_DCT1(nx_global, work_re, -1)
      end select

      phi = 0.0_ISFFT_REAL_KIND
      do i = 1, nx_global
         phi(i,1,1) = work_re(i)
      end do
   end subroutine ISFFT_1D_Mixed_Poisson_serial

   subroutine Algebra_Periodic(spectrum_re, spectrum_im)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: spectrum_re(nx_global), spectrum_im(nx_global)
      type(ISFFT_Scheme_Coefficients) :: coeff
      real(kind=ISFFT_REAL_KIND) :: lambda
      integer :: i, mode_index

      call ISFFT_Get_Scheme_Coefficients(If_scheme, coeff)
      do i = 1, nx_global/2 + 1
         mode_index = i - 1
         lambda = ISFFT_Lambda_For_BC(BC(1), mode_index, nx_global, SLx, If_scheme, coeff)
         if (lambda /= 0.0_ISFFT_REAL_KIND) then
            spectrum_re(i) = -spectrum_re(i)/lambda
            spectrum_im(i) = -spectrum_im(i)/lambda
         else
            spectrum_re(i) = 0.0_ISFFT_REAL_KIND
            spectrum_im(i) = 0.0_ISFFT_REAL_KIND
         end if
      end do
   end subroutine Algebra_Periodic

   subroutine Algebra_Real(spectrum)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: spectrum(nx_global)
      type(ISFFT_Scheme_Coefficients) :: coeff
      real(kind=ISFFT_REAL_KIND) :: lambda
      integer :: i, mode_index

      call ISFFT_Get_Scheme_Coefficients(If_scheme, coeff)
      do i = 1, nx_global
         if (BC(1) == ISFFT_BC_DIRICHLET .and. (i == 1 .or. i == nx_global)) then
            spectrum(i) = 0.0_ISFFT_REAL_KIND
            cycle
         end if
         mode_index = i - 1
         lambda = ISFFT_Lambda_For_BC(BC(1), mode_index, nx_global, SLx, If_scheme, coeff)
         if (lambda /= 0.0_ISFFT_REAL_KIND) then
            spectrum(i) = -spectrum(i)/lambda
         else
            spectrum(i) = 0.0_ISFFT_REAL_KIND
         end if
      end do
   end subroutine Algebra_Real

end module ISFFT_1D_Mixed_Poisson
