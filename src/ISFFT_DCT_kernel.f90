!------------------------------------------------------------------------------
! ISFFT_DCT_kernel
! Scalar and batched DCT-I execution with prepared workspaces.
!------------------------------------------------------------------------------
module ISFFT_DCT_kernel
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   use ISFFT_Workspace, only: fft_precompute_slots, Real_Transform_Workspace, &
      Require_Real_Transform_Precompute, Require_Real_Workspace
   use ISFFT_FFT_kernel, only: ISFFT_1D_FFT_c2c, ISFFT_1D_FFT_c2c_complex_batch, &
      ISFFT_1D_FFT_r2c_complex_batch, ISFFT_1D_FFT_r2c_odd_batch
   implicit none
   private
   public :: ISFFT_1D_DCT1, ISFFT_1D_DCT1_batch, ISFFT_DCT1_direct
   public :: ISFFT_Execute_DCT_Plan

contains

   !------------------------------------------------------------------------------
   !! Scalar and batched cosine transforms
   !------------------------------------------------------------------------------
   !! Here n includes both endpoints. Preparation uses the interval count n-1.
   !! 1D DCT-I interface.
   !
   !! sign >= 0 : forward  DCT-I
   !! sign <  0 : inverse  DCT-I, i.e. x = (1/(2*(n-1))) * DCT-I(y)
   !
   !! Even interval counts use a reduced real FFT. Odd counts use paired
   !! real lines in batches, with the extension FFT retained for scalar calls.
   !! ISFFT_DCT1_direct remains available as the O(n^2) reference kernel.
   subroutine ISFFT_1D_DCT1(n, fr, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout), target :: fr(n)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: batch_view(:, :, :)

      real(kind=ISFFT_REAL_KIND) :: scale
      real(kind=ISFFT_REAL_KIND) :: ur, ui, vr, vi
      real(kind=ISFFT_REAL_KIND) :: t1r, t1i, t2r, t2i, wr, wi, qr, qi
      integer :: k, p, fft_len, packed_len, slot_fft

      if (n < 1) stop 'DCT1: n must be positive'
      call Require_Real_Transform_Precompute(max(1, n - 1), slot_fft)

      if (n == 1) return

      if (mod(n - 1, 2) == 0) then
         batch_view(1:n, 1:1, 1:1) => fr
         call DCT1_Reduced_Batch(n, batch_view, sign, 1, 1)
         return
      end if

      fft_len = 2 * (n - 1)
      packed_len = fft_len / 2
      call Require_Real_Transform_Precompute(packed_len, slot_fft)
      associate (packed_r => fft_precompute_slots(slot_fft)%aux_r(1:packed_len), &
         packed_i => fft_precompute_slots(slot_fft)%aux_i(1:packed_len))
         do k = 0, packed_len - 1
            p = 2 * k
            if (p <= packed_len) then
               packed_r(k + 1) = fr(p + 1)
            else
               packed_r(k + 1) = fr(fft_len - p + 1)
            end if
            p = 2 * k + 1
            if (p <= packed_len) then
               packed_i(k + 1) = fr(p + 1)
            else
               packed_i(k + 1) = fr(fft_len - p + 1)
            end if
         end do

         call ISFFT_1D_FFT_c2c(packed_len, packed_r, packed_i, 1)
         scale = 0.5_ISFFT_REAL_KIND * real(fft_len, ISFFT_REAL_KIND)
         if (sign < 0) scale = 0.5_ISFFT_REAL_KIND
         fr(1) = scale * (packed_r(1) + packed_i(1))
         fr(n) = scale * (packed_r(1) - packed_i(1))
         do k = 1, packed_len - 1
            ur = packed_r(k + 1)
            ui = packed_i(k + 1)
            vr = packed_r(packed_len - k + 1)
            vi = packed_i(packed_len - k + 1)
            t1r = 0.5_ISFFT_REAL_KIND * (ur + vr)
            t1i = 0.5_ISFFT_REAL_KIND * (ui - vi)
            t2r = 0.5_ISFFT_REAL_KIND * (ur - vr)
            t2i = 0.5_ISFFT_REAL_KIND * (ui + vi)
            wr = fft_precompute_slots(slot_fft)%real_twc(k + 1)
            wi = -fft_precompute_slots(slot_fft)%real_tws(k + 1)
            qr = wr * t2r - wi * t2i
            qi = wr * t2i + wi * t2r
            fr(k + 1) = scale * (t1r + qi)
         end do
      end associate
   end subroutine ISFFT_1D_DCT1

   subroutine ISFFT_1D_DCT1_batch(n, fr, sign, j_stop, k_stop)
      integer, intent(in) :: n, sign, j_stop, k_stop
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: packed(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: scale, ur, ui, vr, vi, t1r, t1i, t2r, t2i, wr, wi, qr, qi, value_r, value_i
      integer :: j, k, q, fft_len, packed_len, slot_fft
      type(Real_Transform_Workspace), pointer :: work

      if (n < 1) error stop 'DCT1 batch: n must be positive'
      call Require_Real_Transform_Precompute(max(1, n - 1), slot_fft)
      if (n == 1 .or. j_stop <= 0 .or. k_stop <= 0) return
      !! Even interval counts admit an N-point real FFT (N/2 complex),
      !! avoiding the old N-point complex FFT of the doubled extension.
      !! For odd N, pair independent real lines in one complex FFT.
      if (mod(n - 1, 2) == 0 .or. j_stop * k_stop > 1) then
         call DCT1_Reduced_Batch(n, fr, sign, j_stop, k_stop)
         return
      end if
      fft_len = 2 * (n - 1)
      packed_len = n - 1
      call Require_Real_Transform_Precompute(packed_len, slot_fft)
      call Require_Real_Workspace(packed_len, j_stop * k_stop, work)
      packed(1:packed_len, 1:j_stop, 1:k_stop) => work%batch_storage(1:packed_len * j_stop * k_stop)
      do k = 1, k_stop
         do j = 1, j_stop
            do q = 0, packed_len - 1
               value_r = fr(work%dct_pack_real_index(q + 1), j, k)
               value_i = fr(work%dct_pack_imag_index(q + 1), j, k)
               packed(q + 1, j, k) = cmplx(value_r, value_i, kind=ISFFT_REAL_KIND)
            end do
         end do
      end do
      call ISFFT_1D_FFT_c2c_complex_batch(packed_len, packed, +1)
      scale = 0.5_ISFFT_REAL_KIND * real(fft_len, ISFFT_REAL_KIND)
      if (sign < 0) scale = 0.5_ISFFT_REAL_KIND
      do k = 1, k_stop
         do j = 1, j_stop
            fr(1, j, k) = scale * (real(packed(1, j, k), ISFFT_REAL_KIND) + aimag(packed(1, j, k)))
            fr(n, j, k) = scale * (real(packed(1, j, k), ISFFT_REAL_KIND) - aimag(packed(1, j, k)))
            do q = 1, packed_len - 1
               ur = real(packed(q + 1, j, k), ISFFT_REAL_KIND)
               ui = aimag(packed(q + 1, j, k))
               vr = real(packed(packed_len - q + 1, j, k), ISFFT_REAL_KIND)
               vi = aimag(packed(packed_len - q + 1, j, k))
               t1r = 0.5_ISFFT_REAL_KIND * (ur + vr)
               t1i = 0.5_ISFFT_REAL_KIND * (ui - vi)
               t2r = 0.5_ISFFT_REAL_KIND * (ur - vr)
               t2i = 0.5_ISFFT_REAL_KIND * (ui + vi)
               wr = fft_precompute_slots(slot_fft)%real_twc(q + 1)
               wi = -fft_precompute_slots(slot_fft)%real_tws(q + 1)
               qr = wr * t2r - wi * t2i
               qi = wr * t2i + wi * t2r
               fr(q + 1, j, k) = scale * (t1r + qi)
            end do
         end do
      end do
   end subroutine ISFFT_1D_DCT1_batch

   !! Pairing x(j) and x(N-j) separates the even coefficients from
   !! consecutive odd-coefficient differences.  A prefix sum reconstructs
   !! the latter.  Normalizations match the public direct-reference kernels.
   subroutine DCT1_Reduced_Batch(n, fr, sign, j_stop, k_stop)
      integer, intent(in) :: n, sign, j_stop, k_stop
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(:, :, :)
      integer :: l, m, j, k, q, slot, last_pair
      type(Real_Transform_Workspace), pointer :: work
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: dct_real_work(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: dct_real_spectrum(:, :, :), dct_odd_packed(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: a, b, u, v, acc, scale

      l = n - 1
      m = l / 2
      call Require_Real_Transform_Precompute(l, slot)
      call Require_Real_Workspace(l, j_stop * k_stop, work)
      dct_real_work(1:l, 1:j_stop, 1:k_stop) => work%real_work(1:l * j_stop * k_stop)
      dct_real_spectrum(1:m + 1, 1:j_stop, 1:k_stop) => work%spectrum(1:(m + 1) * j_stop * k_stop)
      dct_odd_packed(1:l, 1:(j_stop * k_stop + 1) / 2, 1:1) => work%odd_packed(1:l * ((j_stop * k_stop + 1) / 2))
      last_pair = (l - 1) / 2
      do k = 1, k_stop
         do j = 1, j_stop
            dct_real_work(1, j, k) = 0.5_ISFFT_REAL_KIND * (fr(1, j, k) + fr(n, j, k))
            if (mod(l, 2) == 0) dct_real_work(m + 1, j, k) = fr(m + 1, j, k)
            acc = 0.5_ISFFT_REAL_KIND * (fr(1, j, k) - fr(n, j, k))
            do q = 1, last_pair
               a = fr(q + 1, j, k)
               b = fr(l - q + 1, j, k)
               u = 0.5_ISFFT_REAL_KIND * (a + b)
               v = fft_precompute_slots(slot)%real_tws(q + 1) * (a - b)
               dct_real_work(q + 1, j, k) = u - v
               dct_real_work(l - q + 1, j, k) = u + v
               acc = acc + fft_precompute_slots(slot)%real_twc(q + 1) * (a - b)
            end do
            fr(2, j, k) = acc
         end do
      end do
      !! ISFFT stores interior R2C modes with twice the endpoint weight.
      if (mod(l, 2) == 0) then
         call ISFFT_1D_FFT_r2c_complex_batch(l, dct_real_work, dct_real_spectrum, +1, .true.)
      else
         call ISFFT_1D_FFT_r2c_odd_batch(l, dct_real_work, dct_real_spectrum, dct_odd_packed)
      end if
      scale = 2.0_ISFFT_REAL_KIND
      if (sign < 0) scale = 1.0_ISFFT_REAL_KIND / real(l, ISFFT_REAL_KIND)
      do k = 1, k_stop
         do j = 1, j_stop
            acc = fr(2, j, k)
            fr(1, j, k) = scale * real(dct_real_spectrum(1, j, k), ISFFT_REAL_KIND)
            fr(2, j, k) = scale * acc
            do q = 1, last_pair
               fr(2 * q + 1, j, k) = 0.5_ISFFT_REAL_KIND * scale * real(dct_real_spectrum(q + 1, j, k), ISFFT_REAL_KIND)
               acc = acc - 0.5_ISFFT_REAL_KIND * aimag(dct_real_spectrum(q + 1, j, k))
               fr(2 * q + 2, j, k) = scale * acc
            end do
            if (mod(l, 2) == 0) fr(n, j, k) = scale * real(dct_real_spectrum(m + 1, j, k), ISFFT_REAL_KIND)
         end do
      end do
   end subroutine DCT1_Reduced_Batch

   !------------------------------------------------------------------------------
   !! Direct reference transform
   !------------------------------------------------------------------------------
   subroutine ISFFT_DCT1_direct(n, x, y)
      integer, intent(in) :: n
      real(kind=ISFFT_REAL_KIND), intent(in) :: x(n)
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
            theta = PI * real((i - 1) * (k - 1), kind=ISFFT_REAL_KIND) / real(n - 1, kind=ISFFT_REAL_KIND)
            y(k) = y(k) + 2.0_ISFFT_REAL_KIND * x(i) * cos(theta)
         end do
      end do
   end subroutine ISFFT_DCT1_direct

   !------------------------------------------------------------------------------
   !! Prepared cosine kernel dispatch
   !------------------------------------------------------------------------------
   subroutine ISFFT_Execute_DCT_Plan(n, values, sign, method, width)
      integer, intent(in) :: n, sign, method, width
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: values(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: v(:, :, :)
      integer :: lines, first, last, j, prepared_slot

      call Require_Real_Transform_Precompute(max(1, n - 1), prepared_slot)
      lines = size(values) / n
      if (lines == 0) return
      v(1:n, 1:lines, 1:1) => values
      if (method == 0) then
         do j = 1, lines
            call ISFFT_1D_DCT1(n, v(:, j, 1), sign)
         end do
      else
         do first = 1, lines, max(1, width)
            last = min(lines, first + max(1, width) - 1)
            call ISFFT_1D_DCT1_batch(n, v(:, first:last, :), sign, last - first + 1, 1)
         end do
      end if
   end subroutine ! ISFFT_Execute_DCT_Plan

end module ISFFT_DCT_kernel
