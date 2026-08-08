module ISFFT_FFT_kernel
   use ISFFT_MPI_Constants
   use ISFFT_FFT_Pre, only: Ensure_FFT_Precompute, Ensure_Rader_Precompute
   use ISFFT_Buffer, only: fft_precompute_slots
   implicit none

   private
   public :: ISFFT_1D_FFT_c2c, ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c_reference, &
      ISFFT_1D_FFT_c2c_rader, ISFFT_1D_FFT_c2c_bluestein

contains

   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_core_radix2(len, fr, fi, sign, rev, twc, tws, work_r, work_i)
      integer, intent(in) :: len, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(len), fi(len)
      integer, intent(in) :: rev(:)
      real(kind=ISFFT_REAL_KIND), intent(in) :: twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: work_r(:), work_i(:)

      integer :: i, j, l, m, iistep, kstep
      real(kind=ISFFT_REAL_KIND) :: wr, wi, tr, ti, xn

      do i = 1, len
         j = rev(i)
         work_r(j) = fr(i)
         work_i(j) = fi(i)
      end do
      fr = work_r(1:len)
      fi = work_i(1:len)

      l = 1
      do while (l < len)
         iistep = 2*l
         kstep = len / iistep
         do m = 0, l-1
            wr = twc(1 + m*kstep)
            wi = -real(sign, kind=ISFFT_REAL_KIND) * tws(1 + m*kstep)
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
         xn = 1.0_ISFFT_REAL_KIND / real(len, kind=ISFFT_REAL_KIND)
         do i = 1, len
            fr(i) = fr(i) * xn
            fi(i) = fi(i) * xn
         end do
      end if
   end subroutine ISFFT_1D_FFT_core_radix2
   !----------------------------------------------------------------------
   recursive subroutine ISFFT_1D_FFT_core_mixed(n_total, n_current, input_offset, input_stride, &
      output_offset, level, fr, fi, out_r, out_i, scratch_r, scratch_i, sign, stage_radix, twc, tws)
      integer, intent(in) :: n_total, n_current, input_offset, input_stride
      integer, intent(in) :: output_offset, level, sign
      real(kind=ISFFT_REAL_KIND), intent(in) :: fr(:), fi(:), twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: out_r(:), out_i(:), scratch_r(:), scratch_i(:)
      integer, intent(in) :: stage_radix(:)

      integer :: radix, child_size, child, k, j, input_index, root_index
      real(kind=ISFFT_REAL_KIND) :: wr, wi, tr, ti, sum_r, sum_i

      if (n_current == 1) then
         out_r(output_offset) = fr(input_offset)
         out_i(output_offset) = fi(input_offset)
         return
      end if

      radix = stage_radix(level)
      child_size = n_current / radix
      do child = 0, radix - 1
         call ISFFT_1D_FFT_core_mixed(n_total, child_size, input_offset + child*input_stride, &
            input_stride*radix, output_offset + child*child_size, level + 1, fr, fi, out_r, out_i, &
            scratch_r, scratch_i, sign, stage_radix, twc, tws)
      end do

      scratch_r(output_offset:output_offset+n_current-1) = out_r(output_offset:output_offset+n_current-1)
      scratch_i(output_offset:output_offset+n_current-1) = out_i(output_offset:output_offset+n_current-1)

      select case (radix)
      case (2)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            input_index = input_index + child_size
            root_index = mod(k*(n_total/n_current), n_total)
            wr = twc(root_index + 1)
            wi = -real(sign, ISFFT_REAL_KIND)*tws(root_index + 1)
            tr = wr*scratch_r(input_index) - wi*scratch_i(input_index)
            ti = wr*scratch_i(input_index) + wi*scratch_r(input_index)
            out_r(output_offset + k) = sum_r + tr
            out_i(output_offset + k) = sum_i + ti
         end do
      case (3)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            do j = 1, 2
               input_index = output_offset + j*child_size + mod(k, child_size)
               root_index = mod(j*k*(n_total/n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND)*tws(root_index + 1)
               tr = wr*scratch_r(input_index) - wi*scratch_i(input_index)
               ti = wr*scratch_i(input_index) + wi*scratch_r(input_index)
               sum_r = sum_r + tr
               sum_i = sum_i + ti
            end do
            out_r(output_offset + k) = sum_r
            out_i(output_offset + k) = sum_i
         end do
      case (5)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            do j = 1, 4
               input_index = output_offset + j*child_size + mod(k, child_size)
               root_index = mod(j*k*(n_total/n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND)*tws(root_index + 1)
               tr = wr*scratch_r(input_index) - wi*scratch_i(input_index)
               ti = wr*scratch_i(input_index) + wi*scratch_r(input_index)
               sum_r = sum_r + tr
               sum_i = sum_i + ti
            end do
            out_r(output_offset + k) = sum_r
            out_i(output_offset + k) = sum_i
         end do
      case (7)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            do j = 1, 6
               input_index = output_offset + j*child_size + mod(k, child_size)
               root_index = mod(j*k*(n_total/n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND)*tws(root_index + 1)
               tr = wr*scratch_r(input_index) - wi*scratch_i(input_index)
               ti = wr*scratch_i(input_index) + wi*scratch_r(input_index)
               sum_r = sum_r + tr
               sum_i = sum_i + ti
            end do
            out_r(output_offset + k) = sum_r
            out_i(output_offset + k) = sum_i
         end do
      case default
         do k = 0, n_current - 1
            sum_r = 0.0_ISFFT_REAL_KIND
            sum_i = 0.0_ISFFT_REAL_KIND
            do j = 0, radix - 1
               input_index = output_offset + j*child_size + mod(k, child_size)
               root_index = mod(j*k*(n_total/n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND)*tws(root_index + 1)
               tr = wr*scratch_r(input_index) - wi*scratch_i(input_index)
               ti = wr*scratch_i(input_index) + wi*scratch_r(input_index)
               sum_r = sum_r + tr
               sum_i = sum_i + ti
            end do
            out_r(output_offset + k) = sum_r
            out_i(output_offset + k) = sum_i
         end do
      end select
   end subroutine ISFFT_1D_FFT_core_mixed
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_core_bluestein(n, fr, fi, sign, slot_idx)
      integer, intent(in) :: n, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: i, convolution_len
      real(kind=ISFFT_REAL_KIND) :: chirp_i, tr, ti, kr, ki, scale

      convolution_len = fft_precompute_slots(slot_idx)%bluestein_m
      associate(br => fft_precompute_slots(slot_idx)%bluestein_r, &
         bi => fft_precompute_slots(slot_idx)%bluestein_i, &
         work_r => fft_precompute_slots(slot_idx)%bluestein_work_r, &
         work_i => fft_precompute_slots(slot_idx)%bluestein_work_i, &
         chirp_c => fft_precompute_slots(slot_idx)%chirp_c, &
         chirp_s => fft_precompute_slots(slot_idx)%chirp_s)

      br = 0.0_ISFFT_REAL_KIND
      bi = 0.0_ISFFT_REAL_KIND
      do i = 1, n
         chirp_i = -real(sign, ISFFT_REAL_KIND)*chirp_s(i)
         br(i) = fr(i)*chirp_c(i) - fi(i)*chirp_i
         bi(i) = fr(i)*chirp_i + fi(i)*chirp_c(i)
      end do

      call ISFFT_1D_FFT_core_radix2(convolution_len, br, bi, +1, &
         fft_precompute_slots(slot_idx)%bluestein_rev, &
         fft_precompute_slots(slot_idx)%bluestein_twc, &
         fft_precompute_slots(slot_idx)%bluestein_tws, work_r, work_i)

      do i = 1, convolution_len
         if (sign >= 0) then
            kr = fft_precompute_slots(slot_idx)%kernel_pos_r(i)
            ki = fft_precompute_slots(slot_idx)%kernel_pos_i(i)
         else
            kr = fft_precompute_slots(slot_idx)%kernel_neg_r(i)
            ki = fft_precompute_slots(slot_idx)%kernel_neg_i(i)
         end if
         tr = br(i)*kr - bi(i)*ki
         ti = br(i)*ki + bi(i)*kr
         br(i) = tr
         bi(i) = ti
      end do

      call ISFFT_1D_FFT_core_radix2(convolution_len, br, bi, -1, &
         fft_precompute_slots(slot_idx)%bluestein_rev, &
         fft_precompute_slots(slot_idx)%bluestein_twc, &
         fft_precompute_slots(slot_idx)%bluestein_tws, work_r, work_i)

      scale = 1.0_ISFFT_REAL_KIND
      if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND/real(n, ISFFT_REAL_KIND)
      do i = 1, n
         chirp_i = -real(sign, ISFFT_REAL_KIND)*chirp_s(i)
         tr = br(i)*chirp_c(i) - bi(i)*chirp_i
         ti = br(i)*chirp_i + bi(i)*chirp_c(i)
         fr(i) = scale*tr
         fi(i) = scale*ti
      end do
      end associate
   end subroutine ISFFT_1D_FFT_core_bluestein
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_rader_convolution_transform(q, input_r, input_i, output_r, output_i, &
      scratch_r, scratch_i, sign, slot_idx)
      integer, intent(in) :: q, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(in) :: input_r(q), input_i(q)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: output_r(q), output_i(q)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: scratch_r(q), scratch_i(q)

      output_r = input_r
      output_i = input_i
      if (fft_precompute_slots(slot_idx)%rader_is_radix2) then
         call ISFFT_1D_FFT_core_radix2(q, output_r, output_i, sign, &
            fft_precompute_slots(slot_idx)%rader_rev, &
            fft_precompute_slots(slot_idx)%rader_twc, fft_precompute_slots(slot_idx)%rader_tws, &
            scratch_r, scratch_i)
      else
         call ISFFT_1D_FFT_core_mixed(q, q, 1, 1, 1, 1, input_r, input_i, &
            output_r, output_i, scratch_r, scratch_i, sign, &
            fft_precompute_slots(slot_idx)%rader_stage_radix, &
            fft_precompute_slots(slot_idx)%rader_twc, fft_precompute_slots(slot_idx)%rader_tws)
         if (sign >= 0) then
            output_r = output_r/real(q, ISFFT_REAL_KIND)
            output_i = output_i/real(q, ISFFT_REAL_KIND)
         end if
      end if
   end subroutine ISFFT_1D_FFT_rader_convolution_transform
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_core_rader(n, fr, fi, sign, slot_idx)
      integer, intent(in) :: n, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: i, q, output_index
      real(kind=ISFFT_REAL_KIND) :: x0r, x0i, sum_r, sum_i, tr, ti, kr, ki, scale

      if (.not. fft_precompute_slots(slot_idx)%rader_ready) then
         error stop 'ISFFT_1D_FFT_core_rader: Rader plan is not available'
      end if
      q = n - 1
      x0r = fr(1)
      x0i = fi(1)
      sum_r = sum(fr)
      sum_i = sum(fi)

      associate(data_r => fft_precompute_slots(slot_idx)%rader_data_r, &
         data_i => fft_precompute_slots(slot_idx)%rader_data_i, &
         work_r1 => fft_precompute_slots(slot_idx)%rader_work_r1, &
         work_i1 => fft_precompute_slots(slot_idx)%rader_work_i1, &
         work_r2 => fft_precompute_slots(slot_idx)%rader_work_r2, &
         work_i2 => fft_precompute_slots(slot_idx)%rader_work_i2)

      do i = 1, q
         data_r(i) = fr(fft_precompute_slots(slot_idx)%rader_perm(i)+1)
         data_i(i) = fi(fft_precompute_slots(slot_idx)%rader_perm(i)+1)
      end do
      call ISFFT_1D_FFT_rader_convolution_transform(q, data_r, data_i, work_r1, work_i1, &
         work_r2, work_i2, +1, slot_idx)

      do i = 1, q
         if (sign >= 0) then
            kr = fft_precompute_slots(slot_idx)%rader_kernel_pos_r(i)
            ki = fft_precompute_slots(slot_idx)%rader_kernel_pos_i(i)
         else
            kr = fft_precompute_slots(slot_idx)%rader_kernel_neg_r(i)
            ki = fft_precompute_slots(slot_idx)%rader_kernel_neg_i(i)
         end if
         tr = work_r1(i)*kr - work_i1(i)*ki
         ti = work_r1(i)*ki + work_i1(i)*kr
         work_r1(i) = tr
         work_i1(i) = ti
      end do

      call ISFFT_1D_FFT_rader_convolution_transform(q, work_r1, work_i1, data_r, data_i, &
         work_r2, work_i2, -1, slot_idx)

      scale = 1.0_ISFFT_REAL_KIND
      if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND/real(n, ISFFT_REAL_KIND)
      fr(1) = scale*sum_r
      fi(1) = scale*sum_i
      do i = 1, q
         output_index = fft_precompute_slots(slot_idx)%rader_out_perm(i) + 1
         fr(output_index) = scale*(x0r + data_r(i))
         fi(output_index) = scale*(x0i + data_i(i))
      end do
      end associate
   end subroutine ISFFT_1D_FFT_core_rader
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_execute_reference(n, fr, fi, sign, slot_idx)
      integer, intent(in) :: n, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      real(kind=ISFFT_REAL_KIND) :: scale

      if (fft_precompute_slots(slot_idx)%is_radix2) then
         call ISFFT_1D_FFT_core_radix2(n, fr, fi, sign, fft_precompute_slots(slot_idx)%rev, &
            fft_precompute_slots(slot_idx)%twc, fft_precompute_slots(slot_idx)%tws, &
            fft_precompute_slots(slot_idx)%work_r1, fft_precompute_slots(slot_idx)%work_i1)
         return
      end if

      call ISFFT_1D_FFT_core_mixed(n, n, 1, 1, 1, 1, fr, fi, &
         fft_precompute_slots(slot_idx)%work_r1, fft_precompute_slots(slot_idx)%work_i1, &
         fft_precompute_slots(slot_idx)%work_r2, fft_precompute_slots(slot_idx)%work_i2, sign, &
         fft_precompute_slots(slot_idx)%stage_radix, fft_precompute_slots(slot_idx)%twc, &
         fft_precompute_slots(slot_idx)%tws)
      fr = fft_precompute_slots(slot_idx)%work_r1
      fi = fft_precompute_slots(slot_idx)%work_i1

      if (sign >= 0) then
         scale = 1.0_ISFFT_REAL_KIND/real(n, ISFFT_REAL_KIND)
         fr = scale*fr
         fi = scale*fi
      end if
   end subroutine ISFFT_1D_FFT_execute_reference
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_execute(n, fr, fi, sign, slot_idx)
      integer, intent(in) :: n, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)

      if (fft_precompute_slots(slot_idx)%is_radix2) then
         call ISFFT_1D_FFT_core_radix2(n, fr, fi, sign, fft_precompute_slots(slot_idx)%rev, &
            fft_precompute_slots(slot_idx)%twc, fft_precompute_slots(slot_idx)%tws, &
            fft_precompute_slots(slot_idx)%work_r1, fft_precompute_slots(slot_idx)%work_i1)
         return
      end if

      if (fft_precompute_slots(slot_idx)%use_bluestein) then
         call ISFFT_1D_FFT_core_bluestein(n, fr, fi, sign, slot_idx)
      else if (fft_precompute_slots(slot_idx)%use_rader) then
         call ISFFT_1D_FFT_core_rader(n, fr, fi, sign, slot_idx)
      else
         call ISFFT_1D_FFT_execute_reference(n, fr, fi, sign, slot_idx)
      end if
   end subroutine ISFFT_1D_FFT_execute
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Ensure_FFT_Precompute(n, slot_n)
      call ISFFT_1D_FFT_execute(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c_reference(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Ensure_FFT_Precompute(n, slot_n)
      call ISFFT_1D_FFT_execute_reference(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_reference
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c_rader(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Ensure_Rader_Precompute(n, slot_n)
      call ISFFT_1D_FFT_core_rader(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_rader
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c_bluestein(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Ensure_FFT_Precompute(n, slot_n)
      call ISFFT_1D_FFT_core_bluestein(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_bluestein
   !----------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_r2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)

      integer :: M, k, slot_n, slot_m
      real(kind=ISFFT_REAL_KIND) :: ur, ui, vr, vi
      real(kind=ISFFT_REAL_KIND) :: t1r, t1i, t2r, t2i
      real(kind=ISFFT_REAL_KIND) :: wr, wi, qr, qi
      real(kind=ISFFT_REAL_KIND) :: a0, b0

      if (mod(n,2) /= 0) stop 'FFT1D_r2c: n must be even'
      M = n/2

      call Ensure_FFT_Precompute(n, slot_n)
      call Ensure_FFT_Precompute(M, slot_m)

      associate(yr => fft_precompute_slots(slot_n)%work_r1(1:M), &
         yi => fft_precompute_slots(slot_n)%work_i1(1:M))
      if (sign >= 0) then
         do k = 1, M
            yr(k) = fr(2*k - 1)
            yi(k) = fr(2*k)
         end do

         call ISFFT_1D_FFT_execute(M, yr, yi, +1, slot_m)

         a0 = yr(1)
         b0 = yi(1)
         fr(1)   = 0.5_ISFFT_REAL_KIND*(a0 + b0)
         fi(1)   = 0.0_ISFFT_REAL_KIND
         fr(M+1) = 0.5_ISFFT_REAL_KIND*(a0 - b0)
         fi(M+1) = 0.0_ISFFT_REAL_KIND

         do k = 1, M-1
            ur = yr(k+1)
            ui = yi(k+1)
            vr = yr(M-k+1)
            vi = yi(M-k+1)

            t1r = 0.5_ISFFT_REAL_KIND * (ur + vr)
            t1i = 0.5_ISFFT_REAL_KIND * (ui - vi)
            t2r = 0.5_ISFFT_REAL_KIND * (ur - vr)
            t2i = 0.5_ISFFT_REAL_KIND * (ui + vi)

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
            t1r = 0.5_ISFFT_REAL_KIND*(fr(k+1) + fr(M-k+1))
            t1i = 0.5_ISFFT_REAL_KIND*(fi(k+1) - fi(M-k+1))
            t2r = 0.5_ISFFT_REAL_KIND*(fr(k+1) - fr(M-k+1))
            t2i = 0.5_ISFFT_REAL_KIND*(fi(k+1) + fi(M-k+1))

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

         call ISFFT_1D_FFT_execute(M, yr, yi, -1, slot_m)

         do k = 1, M
            fr(2*k - 1) = yr(k)
            fr(2*k) = yi(k)
         end do
      end if
      end associate

   end subroutine ISFFT_1D_FFT_r2c
   !----------------------------------------------------------------------

end module ISFFT_FFT_kernel
