module ISFFT_DCT_kernel
   use ISFFT_MPI_Constants
   use ISFFT_FFT_Pre, only: Ensure_Real_Transform_Precompute
   use ISFFT_FFT_kernel, only: ISFFT_1D_FFT_c2c
   use ISFFT_Buffer, only: fft_precompute_slots
   implicit none

   private
   public :: ISFFT_1D_DCT1, ISFFT_DCT1_direct

contains

   !----------------------------------------------------------------------
   ! 1D DCT-I interface.
   !
   ! sign >= 0 : forward  DCT-I
   ! sign <  0 : inverse  DCT-I, i.e. x = (1/(2*(n-1))) * DCT-I(y)
   !
   ! The production path uses an even extension and a general-length FFT.
   ! ISFFT_DCT1_direct remains available as the O(n^2) reference kernel.
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_DCT1(n, fr, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n)

      real(kind=ISFFT_REAL_KIND) :: scale
      real(kind=ISFFT_REAL_KIND) :: ur, ui, vr, vi
      real(kind=ISFFT_REAL_KIND) :: t1r, t1i, t2r, t2i, wr, wi, qr, qi
      integer :: k, p, fft_len, packed_len, slot_fft

      if (n < 1) stop 'DCT1: n must be positive'

      if (n == 1) return

      fft_len = 2*(n - 1)
      packed_len = fft_len/2
      call Ensure_Real_Transform_Precompute(packed_len, slot_fft)
      associate(packed_r => fft_precompute_slots(slot_fft)%aux_r(1:packed_len), &
         packed_i => fft_precompute_slots(slot_fft)%aux_i(1:packed_len))
      do k = 0, packed_len - 1
         p = 2*k
         if (p <= packed_len) then
            packed_r(k+1) = fr(p+1)
         else
            packed_r(k+1) = fr(fft_len-p+1)
         end if
         p = 2*k + 1
         if (p <= packed_len) then
            packed_i(k+1) = fr(p+1)
         else
            packed_i(k+1) = fr(fft_len-p+1)
         end if
      end do

      call ISFFT_1D_FFT_c2c(packed_len, packed_r, packed_i, 1)
      scale = 0.5_ISFFT_REAL_KIND*real(fft_len, ISFFT_REAL_KIND)
      fr(1) = scale*(packed_r(1) + packed_i(1))
      fr(n) = scale*(packed_r(1) - packed_i(1))
      do k = 1, packed_len - 1
         ur = packed_r(k+1)
         ui = packed_i(k+1)
         vr = packed_r(packed_len-k+1)
         vi = packed_i(packed_len-k+1)
         t1r = 0.5_ISFFT_REAL_KIND*(ur + vr)
         t1i = 0.5_ISFFT_REAL_KIND*(ui - vi)
         t2r = 0.5_ISFFT_REAL_KIND*(ur - vr)
         t2i = 0.5_ISFFT_REAL_KIND*(ui + vi)
         wr = fft_precompute_slots(slot_fft)%real_twc(k+1)
         wi = -fft_precompute_slots(slot_fft)%real_tws(k+1)
         qr = wr*t2r - wi*t2i
         qi = wr*t2i + wi*t2r
         fr(k+1) = scale*(t1r + qi)
      end do
      if (sign < 0) fr = fr/real(fft_len, ISFFT_REAL_KIND)
      end associate
   end subroutine ISFFT_1D_DCT1
   !----------------------------------------------------------------------
   subroutine ISFFT_DCT1_direct(n, x, y)
      integer, intent(in) :: n
      real(kind=ISFFT_REAL_KIND), intent(in)  :: x(n)
      real(kind=ISFFT_REAL_KIND), intent(out) :: y(n)

      integer :: i, k
      real(kind=ISFFT_REAL_KIND) :: theta, endpoint_sign
      real(kind=ISFFT_REAL_KIND), parameter :: PI = acos(-1.0_ISFFT_REAL_KIND)

      y = 0.0_ISFFT_REAL_KIND
      do k = 1, n
         if (mod(k - 1, 2) == 0) then
            endpoint_sign = 1.0_ISFFT_REAL_KIND
         else
            endpoint_sign = -1.0_ISFFT_REAL_KIND
         end if

         y(k) = x(1) + endpoint_sign * x(n)
         do i = 2, n - 1
            theta = PI * real((i - 1)*(k - 1), kind=ISFFT_REAL_KIND) &
               / real(n - 1, kind=ISFFT_REAL_KIND)
            y(k) = y(k) + 2.0_ISFFT_REAL_KIND * x(i) * cos(theta)
         end do
      end do
   end subroutine ISFFT_DCT1_direct
   !----------------------------------------------------------------------

end module ISFFT_DCT_kernel
