module SPLASH_FFT_kernel
   use SPLASH_MPI_Constants
   use SPLASH_FFT_Pre, only: Ensure_FFT_Precompute
   use SPLASH_Buffer, only: fft_precompute_slots
   implicit none

   private
   public :: SPLASH_1D_FFT_c2c, SPLASH_1D_FFT_r2c

contains

   !----------------------------------------------------------------------
   subroutine SPLASH_1D_FFT_core_radix2(len, fr, fi, sign, rev, twc, tws)
      integer, intent(in) :: len, sign
      real(kind=SPLASH_REAL_KIND), intent(inout) :: fr(len), fi(len)
      integer, intent(in) :: rev(:)
      real(kind=SPLASH_REAL_KIND), intent(in) :: twc(:), tws(:)

      integer :: i, j, l, m, iistep, kstep
      real(kind=SPLASH_REAL_KIND) :: wr, wi, tr, ti, xn
      real(kind=SPLASH_REAL_KIND), allocatable :: rfr(:), rfi(:)

      allocate(rfr(len), rfi(len))
      do i = 1, len
         j = rev(i)
         rfr(j) = fr(i)
         rfi(j) = fi(i)
      end do
      fr = rfr
      fi = rfi
      deallocate(rfr, rfi)

      l = 1
      do while (l < len)
         iistep = 2*l
         kstep = len / iistep
         do m = 0, l-1
            wr = twc(1 + m*kstep)
            wi = -real(sign, kind=SPLASH_REAL_KIND) * tws(1 + m*kstep)
            do i = m+1, len, iistep
               j  = i + l
               tr = wr*fr(j) - wi*fi(j)
               ti = wr*fi(j) + wi*fr(j)
               fr(j) = fr(i) - tr
               fi(j) = fi(i) - ti
               fr(i) = fr(i) + tr
               fi(i) = fi(i) + ti
            end do
         end do
         l = iistep
      end do

      if (sign >= 0) then
         xn = 1.0_SPLASH_REAL_KIND / real(len, kind=SPLASH_REAL_KIND)
         do i = 1, len
            fr(i) = fr(i) * xn
            fi(i) = fi(i) * xn
         end do
      end if
   end subroutine SPLASH_1D_FFT_core_radix2
   !----------------------------------------------------------------------
   subroutine SPLASH_1D_FFT_c2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=SPLASH_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Ensure_FFT_Precompute(n, slot_n)
      call SPLASH_1D_FFT_core_radix2(n, fr, fi, sign, fft_precompute_slots(slot_n)%rev, &
         fft_precompute_slots(slot_n)%twc, fft_precompute_slots(slot_n)%tws)
   end subroutine SPLASH_1D_FFT_c2c
   !----------------------------------------------------------------------
   subroutine SPLASH_1D_FFT_r2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=SPLASH_REAL_KIND), intent(inout) :: fr(n), fi(n)

      integer :: M, k, slot_n, slot_m
      real(kind=SPLASH_REAL_KIND), allocatable :: yr(:), yi(:)
      real(kind=SPLASH_REAL_KIND) :: ur, ui, vr, vi
      real(kind=SPLASH_REAL_KIND) :: t1r, t1i, t2r, t2i
      real(kind=SPLASH_REAL_KIND) :: wr, wi, qr, qi
      real(kind=SPLASH_REAL_KIND) :: a0, b0

      if (mod(n,2) /= 0) stop 'FFT1D_r2c: n must be even'
      M = n/2

      call Ensure_FFT_Precompute(n, slot_n)
      call Ensure_FFT_Precompute(M, slot_m)

      allocate(yr(M), yi(M))

      if (sign >= 0) then
         do k = 1, M
            yr(k) = fr(2*k - 1)
            yi(k) = fr(2*k)
         end do

         call SPLASH_1D_FFT_core_radix2(M, yr, yi, +1, fft_precompute_slots(slot_m)%rev, &
            fft_precompute_slots(slot_m)%twc, fft_precompute_slots(slot_m)%tws)

         a0 = yr(1)
         b0 = yi(1)
         fr(1)   = 0.5_SPLASH_REAL_KIND*(a0 + b0)
         fi(1)   = 0.0_SPLASH_REAL_KIND
         fr(M+1) = 0.5_SPLASH_REAL_KIND*(a0 - b0)
         fi(M+1) = 0.0_SPLASH_REAL_KIND

         do k = 1, M-1
            ur = yr(k+1)
            ui = yi(k+1)
            vr = yr(M-k+1)
            vi = yi(M-k+1)

            t1r = 0.5_SPLASH_REAL_KIND * (ur + vr)
            t1i = 0.5_SPLASH_REAL_KIND * (ui - vi)
            t2r = 0.5_SPLASH_REAL_KIND * (ur - vr)
            t2i = 0.5_SPLASH_REAL_KIND * (ui + vi)

            wr = fft_precompute_slots(slot_n)%twc(k+1)
            wi = -fft_precompute_slots(slot_n)%tws(k+1)

            qr = wr*t2r - wi*t2i
            qi = wr*t2i + wi*t2r

            fr(k+1) = t1r + qi
            fi(k+1) = t1i - qr
         end do

      else
         a0 = fr(1)
         b0 = fr(M+1)
         yr(1) = a0 + b0
         yi(1) = a0 - b0

         do k = 1, M-1
            t1r = 0.5_SPLASH_REAL_KIND*(fr(k+1) + fr(M-k+1))
            t1i = 0.5_SPLASH_REAL_KIND*(fi(k+1) - fi(M-k+1))
            t2r = 0.5_SPLASH_REAL_KIND*(fr(k+1) - fr(M-k+1))
            t2i = 0.5_SPLASH_REAL_KIND*(fi(k+1) + fi(M-k+1))

            wr = fft_precompute_slots(slot_n)%twc(k+1)
            wi = fft_precompute_slots(slot_n)%tws(k+1)

            qr = wr*t2r - wi*t2i
            qi = wr*t2i + wi*t2r

            t2r = -qi
            t2i = qr

            yr(k+1)   = t1r + t2r
            yi(k+1)   = t1i + t2i
            yr(M-k+1) = t1r - t2r
            yi(M-k+1) = -(t1i - t2i)
         end do

         call SPLASH_1D_FFT_core_radix2(M, yr, yi, -1, fft_precompute_slots(slot_m)%rev, &
            fft_precompute_slots(slot_m)%twc, fft_precompute_slots(slot_m)%tws)

         do k = 1, M
            fr(2*k - 1) = yr(k)
            fr(2*k) = yi(k)
         end do
      end if

      deallocate(yr, yi)
   end subroutine SPLASH_1D_FFT_r2c
   !----------------------------------------------------------------------

end module SPLASH_FFT_kernel
