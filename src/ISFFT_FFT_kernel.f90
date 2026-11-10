!------------------------------------------------------------------------------
! ISFFT_FFT_kernel
! Prepared scalar, batched and strided FFT execution kernels.
!------------------------------------------------------------------------------
module ISFFT_FFT_kernel
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND
   use ISFFT_Workspace, only: fft_precompute_slots, Require_FFT_Precompute, Require_Rader_Precompute
   implicit none
   private
   public :: ISFFT_1D_FFT_c2c, ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c_reference, ISFFT_1D_FFT_c2c_rader, &
      ISFFT_1D_FFT_c2c_bluestein, ISFFT_1D_FFT_c2c_batch
   public :: ISFFT_1D_FFT_c2c_complex_batch
   public :: ISFFT_1D_FFT_c2c_rader_batch
   public :: ISFFT_Execute_Axis_Lines
   public :: ISFFT_1D_FFT_c2c_axis_batch
   public :: ISFFT_1D_FFT_r2c_complex_batch
   public :: ISFFT_1D_FFT_r2c_odd_batch
   public :: ISFFT_R2C_Selected_Deferred_Pair_Scale
   public :: ISFFT_Execute_Complex_Plan, ISFFT_Execute_R2C_Plan, ISFFT_Execute_R2C_Forward
   public :: ISFFT_FFT_Validation_Bound
   public :: ISFFT_Execute_FFT_Complex_Strategy, ISFFT_1D_FFT_execute_reference, ISFFT_Execute_R2C_Scalar, &
      ISFFT_Execute_R2C_Batch, ISFFT_Execute_R2C_Generic_Batch, ISFFT_Execute_R2C_Odd_Batch, &
      ISFFT_Is_Three_Times_Power_Of_Two

contains

   !------------------------------------------------------------------------------
   !! Prepared execution entry points
   !------------------------------------------------------------------------------
   subroutine ISFFT_Execute_R2C_Forward(n, source, spectrum, method, width)
      integer, intent(in) :: n, method, width
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: source(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: unused(0, 0, 0)

      select case (method)
      case (5)
         call ISFFT_Execute_R2C_Odd_Batch(n, unused, spectrum, 1, .false., width, source)
      case (2, 4)
         call ISFFT_Execute_R2C_Batch(n, unused, spectrum, 1, width, source=source)
      case (3)
         call ISFFT_Execute_R2C_Generic_Batch(n, unused, spectrum, 1, width, source=source)
      case default
         call ISFFT_Execute_R2C_Scalar(n, unused, spectrum, 1, source)
      end select
   end subroutine ! ISFFT_Execute_R2C_Forward

   !! Method IDs are shared with Transform_Preparation candidate selection.
   !! Update both dispatch and preparation when adding or renumbering methods.
   !! C2C methods: scalar native=1, scalar radix4=2, batch radix4=3,
   !! batch radix2=4, batch Bluestein=5, scalar mixed=6, batch mixed=7/9,
   !! batch Rader=8, batch radix3x2=10, scalar Bluestein=11, scalar Rader=12.
   !! R2C: scalar=1, packed batch=2/4, generic batch=3, odd paired batch=5.
   !! Immutable strategy entry points: all measurements belong to Planner.
   subroutine ISFFT_Execute_Complex_Plan(n, values, sign, method, width)
      integer, intent(in) :: n, sign, method, width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: values(:, :, :)
      integer :: slot

      call Require_FFT_Precompute(n, slot)
      call ISFFT_Execute_FFT_Complex_Strategy(n, values, sign, slot, method, width)
   end subroutine ! ISFFT_Execute_Complex_Plan

   subroutine ISFFT_Execute_R2C_Plan(n, re, spectrum, sign, method, width)
      integer, intent(in) :: n, sign, method, width
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)

      select case (method)
      case (5)
         call ISFFT_Execute_R2C_Odd_Batch(n, re, spectrum, sign, .false., width)
      case (2, 4)
         call ISFFT_Execute_R2C_Batch(n, re, spectrum, sign, width)
      case (3)
         call ISFFT_Execute_R2C_Generic_Batch(n, re, spectrum, sign, width)
      case default
         call ISFFT_Execute_R2C_Scalar(n, re, spectrum, sign)
      end select
   end subroutine ! ISFFT_Execute_R2C_Plan

   !------------------------------------------------------------------------------
   !! Validation bounds and kernel predicates
   !------------------------------------------------------------------------------
   pure real(kind=ISFFT_REAL_KIND) function ISFFT_FFT_Validation_Bound(n) result(bound)
      integer, intent(in) :: n
      real(kind=ISFFT_REAL_KIND) :: accumulated_roundoff

      !! A complex multiply/add needs at most eight real rounded operations.
      !! gamma_k = k*u/(1-k*u) is the standard first-order accumulation bound.
      !! Comparing two independently rounded transforms doubles that bound.
      accumulated_roundoff = 8.0_ISFFT_REAL_KIND * real(max(1, n), ISFFT_REAL_KIND) * epsilon(1.0_ISFFT_REAL_KIND)
      if (accumulated_roundoff >= 1.0_ISFFT_REAL_KIND) then
         bound = huge(bound)
      else
         bound = 2.0_ISFFT_REAL_KIND * accumulated_roundoff / (1.0_ISFFT_REAL_KIND - accumulated_roundoff)
      end if
   end function ISFFT_FFT_Validation_Bound

   logical function ISFFT_R2C_Selected_Deferred_Pair_Scale(n, total_lines) result(selected)
      integer, intent(in) :: n, total_lines
      integer :: slot_n, entry, method

      call Require_FFT_Precompute(n, slot_n)
      method = 1
      do entry = 1, fft_precompute_slots(slot_n)%r2c_strategy_count
         if (fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) == total_lines) then
            method = fft_precompute_slots(slot_n)%r2c_strategy_method(entry)
            exit
         end if
      end do
      selected = method == 2 .or. method == 4
   end function ISFFT_R2C_Selected_Deferred_Pair_Scale

   logical function ISFFT_Is_Three_Times_Power_Of_Two(n) result(matches)
      integer, intent(in) :: n
      integer :: value

      if (n < 3 .or. mod(n, 3) /= 0) then
         matches = .false.
         return
      end if
      value = n / 3
      do while (value > 1 .and. mod(value, 2) == 0)
         value = value / 2
      end do
      matches = value == 1
   end function ISFFT_Is_Three_Times_Power_Of_Two

   !------------------------------------------------------------------------------
   !! Split-real scalar FFT algorithms
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_core_radix2(len, fr, fi, sign, rev, twc, tws, work_r, work_i)
      integer, intent(in) :: len, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(len), fi(len)
      integer, intent(in) :: rev(:)
      real(kind=ISFFT_REAL_KIND), intent(in) :: twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: work_r(:), work_i(:)

      integer :: i, ii, j, l, m, iistep, kstep
      real(kind=ISFFT_REAL_KIND) :: wr, wi, tr, ti, xn

      if (sign >= 0) then
         xn = 1.0_ISFFT_REAL_KIND / real(len, kind=ISFFT_REAL_KIND)
         do i = 1, len
            j = rev(i)
            work_r(j) = xn * fr(i)
            work_i(j) = xn * fi(i)
         end do
      else
         do i = 1, len
            j = rev(i)
            work_r(j) = fr(i)
            work_i(j) = fi(i)
         end do
      end if
      fr = work_r(1:len)
      fi = work_i(1:len)

      l = 1
      do while (l < len)
         iistep = 2 * l
         kstep = len / iistep
         do i = 1, len, iistep
            do m = 0, l - 1
               ii = i + m
               j = ii + l
               wr = twc(1 + m * kstep)
               wi = -real(sign, kind=ISFFT_REAL_KIND) * tws(1 + m * kstep)
               tr = wr * fr(j) - wi * fi(j)
               ti = wr * fi(j) + wi * fr(j)
               fr(j) = fr(ii) - tr
               fi(j) = fi(ii) - ti
               fr(ii) = fr(ii) + tr
               fi(ii) = fi(ii) + ti
            end do
         end do
         l = iistep
      end do

   end subroutine ISFFT_1D_FFT_core_radix2

   recursive subroutine ISFFT_1D_FFT_core_mixed( &
      n_total, n_current, input_offset, input_stride, output_offset, level, &
      fr, fi, out_r, out_i, scratch_r, scratch_i, &
      sign, stage_radix, twc, tws)
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
         call ISFFT_1D_FFT_core_mixed( &
            n_total, child_size, input_offset + child * input_stride, input_stride * radix, &
            output_offset + child * child_size, level + 1, &
            fr, fi, out_r, out_i, scratch_r, scratch_i, &
            sign, stage_radix, twc, tws)
      end do

      scratch_r(output_offset:output_offset + n_current - 1) = out_r(output_offset:output_offset + n_current - 1)
      scratch_i(output_offset:output_offset + n_current - 1) = out_i(output_offset:output_offset + n_current - 1)

      select case (radix)
      case (2)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            input_index = input_index + child_size
            root_index = mod(k * (n_total / n_current), n_total)
            wr = twc(root_index + 1)
            wi = -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
            tr = wr * scratch_r(input_index) - wi * scratch_i(input_index)
            ti = wr * scratch_i(input_index) + wi * scratch_r(input_index)
            out_r(output_offset + k) = sum_r + tr
            out_i(output_offset + k) = sum_i + ti
         end do
      case (3)
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            sum_r = scratch_r(input_index)
            sum_i = scratch_i(input_index)
            do j = 1, 2
               input_index = output_offset + j * child_size + mod(k, child_size)
               root_index = mod(j * k * (n_total / n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
               tr = wr * scratch_r(input_index) - wi * scratch_i(input_index)
               ti = wr * scratch_i(input_index) + wi * scratch_r(input_index)
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
               input_index = output_offset + j * child_size + mod(k, child_size)
               root_index = mod(j * k * (n_total / n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
               tr = wr * scratch_r(input_index) - wi * scratch_i(input_index)
               ti = wr * scratch_i(input_index) + wi * scratch_r(input_index)
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
               input_index = output_offset + j * child_size + mod(k, child_size)
               root_index = mod(j * k * (n_total / n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
               tr = wr * scratch_r(input_index) - wi * scratch_i(input_index)
               ti = wr * scratch_i(input_index) + wi * scratch_r(input_index)
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
               input_index = output_offset + j * child_size + mod(k, child_size)
               root_index = mod(j * k * (n_total / n_current), n_total)
               wr = twc(root_index + 1)
               wi = -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
               tr = wr * scratch_r(input_index) - wi * scratch_i(input_index)
               ti = wr * scratch_i(input_index) + wi * scratch_r(input_index)
               sum_r = sum_r + tr
               sum_i = sum_i + ti
            end do
            out_r(output_offset + k) = sum_r
            out_i(output_offset + k) = sum_i
         end do
      end select
   end subroutine ISFFT_1D_FFT_core_mixed

   subroutine ISFFT_1D_FFT_core_bluestein(n, fr, fi, sign, slot_idx)
      integer, intent(in) :: n, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: i, convolution_len
      real(kind=ISFFT_REAL_KIND) :: chirp_i, tr, ti, kr, ki, scale

      convolution_len = fft_precompute_slots(slot_idx)%bluestein_m
      associate (br => fft_precompute_slots(slot_idx)%bluestein_r, bi => fft_precompute_slots(slot_idx)%bluestein_i, &
         work_r => fft_precompute_slots(slot_idx)%bluestein_work_r, &
         work_i => fft_precompute_slots(slot_idx)%bluestein_work_i, &
         chirp_c => fft_precompute_slots(slot_idx)%chirp_c, chirp_s => fft_precompute_slots(slot_idx)%chirp_s)

         br = 0.0_ISFFT_REAL_KIND
         bi = 0.0_ISFFT_REAL_KIND
         do i = 1, n
            chirp_i = -real(sign, ISFFT_REAL_KIND) * chirp_s(i)
            br(i) = fr(i) * chirp_c(i) - fi(i) * chirp_i
            bi(i) = fr(i) * chirp_i + fi(i) * chirp_c(i)
         end do

         call ISFFT_1D_FFT_core_radix2(convolution_len, br, bi, +1, fft_precompute_slots(slot_idx)%bluestein_rev, &
            fft_precompute_slots(slot_idx)%bluestein_twc, fft_precompute_slots(slot_idx)%bluestein_tws, work_r, work_i)

         do i = 1, convolution_len
            if (sign >= 0) then
               kr = fft_precompute_slots(slot_idx)%kernel_pos_r(i)
               ki = fft_precompute_slots(slot_idx)%kernel_pos_i(i)
            else
               kr = fft_precompute_slots(slot_idx)%kernel_neg_r(i)
               ki = fft_precompute_slots(slot_idx)%kernel_neg_i(i)
            end if
            tr = br(i) * kr - bi(i) * ki
            ti = br(i) * ki + bi(i) * kr
            br(i) = tr
            bi(i) = ti
         end do

         call ISFFT_1D_FFT_core_radix2(convolution_len, br, bi, -1, fft_precompute_slots(slot_idx)%bluestein_rev, &
            fft_precompute_slots(slot_idx)%bluestein_twc, fft_precompute_slots(slot_idx)%bluestein_tws, work_r, work_i)

         scale = 1.0_ISFFT_REAL_KIND
         if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
         do i = 1, n
            chirp_i = -real(sign, ISFFT_REAL_KIND) * chirp_s(i)
            tr = br(i) * chirp_c(i) - bi(i) * chirp_i
            ti = br(i) * chirp_i + bi(i) * chirp_c(i)
            fr(i) = scale * tr
            fi(i) = scale * ti
         end do
      end associate
   end subroutine ISFFT_1D_FFT_core_bluestein

   subroutine ISFFT_1D_FFT_rader_convolution_transform( &
      q, input_r, input_i, &
      output_r, output_i, scratch_r, scratch_i, &
      sign, slot_idx)
      integer, intent(in) :: q, sign, slot_idx
      real(kind=ISFFT_REAL_KIND), intent(in) :: input_r(q), input_i(q)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: output_r(q), output_i(q)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: scratch_r(q), scratch_i(q)

      output_r = input_r
      output_i = input_i
      if (fft_precompute_slots(slot_idx)%rader_is_radix2) then
         call ISFFT_1D_FFT_core_radix2(q, output_r, output_i, sign, fft_precompute_slots(slot_idx)%rader_rev, &
            fft_precompute_slots(slot_idx)%rader_twc, fft_precompute_slots(slot_idx)%rader_tws, scratch_r, scratch_i)
      else
         call ISFFT_1D_FFT_core_mixed( &
            q, q, 1, 1, 1, 1, &
            input_r, input_i, output_r, output_i, scratch_r, scratch_i, &
            sign, fft_precompute_slots(slot_idx)%rader_stage_radix, fft_precompute_slots(slot_idx)%rader_twc, &
            fft_precompute_slots(slot_idx)%rader_tws)
         if (sign >= 0) then
            output_r = output_r / real(q, ISFFT_REAL_KIND)
            output_i = output_i / real(q, ISFFT_REAL_KIND)
         end if
      end if
   end subroutine ISFFT_1D_FFT_rader_convolution_transform

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

      associate (data_r => fft_precompute_slots(slot_idx)%rader_data_r, &
         data_i => fft_precompute_slots(slot_idx)%rader_data_i, &
         work_r1 => fft_precompute_slots(slot_idx)%rader_work_r1, &
         work_i1 => fft_precompute_slots(slot_idx)%rader_work_i1, &
         work_r2 => fft_precompute_slots(slot_idx)%rader_work_r2, &
         work_i2 => fft_precompute_slots(slot_idx)%rader_work_i2)

         do i = 1, q
            data_r(i) = fr(fft_precompute_slots(slot_idx)%rader_perm(i) + 1)
            data_i(i) = fi(fft_precompute_slots(slot_idx)%rader_perm(i) + 1)
         end do
         call ISFFT_1D_FFT_rader_convolution_transform( &
            q, data_r, data_i, &
            work_r1, work_i1, work_r2, work_i2, &
            +1, slot_idx)

         do i = 1, q
            if (sign >= 0) then
               kr = fft_precompute_slots(slot_idx)%rader_kernel_pos_r(i)
               ki = fft_precompute_slots(slot_idx)%rader_kernel_pos_i(i)
            else
               kr = fft_precompute_slots(slot_idx)%rader_kernel_neg_r(i)
               ki = fft_precompute_slots(slot_idx)%rader_kernel_neg_i(i)
            end if
            tr = work_r1(i) * kr - work_i1(i) * ki
            ti = work_r1(i) * ki + work_i1(i) * kr
            work_r1(i) = tr
            work_i1(i) = ti
         end do

         call ISFFT_1D_FFT_rader_convolution_transform( &
            q, work_r1, work_i1, &
            data_r, data_i, work_r2, work_i2, &
            -1, slot_idx)

         scale = 1.0_ISFFT_REAL_KIND
         if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
         fr(1) = scale * sum_r
         fi(1) = scale * sum_i
         do i = 1, q
            output_index = fft_precompute_slots(slot_idx)%rader_out_perm(i) + 1
            fr(output_index) = scale * (x0r + data_r(i))
            fi(output_index) = scale * (x0i + data_i(i))
         end do
      end associate
   end subroutine ISFFT_1D_FFT_core_rader

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

      call ISFFT_1D_FFT_core_mixed( &
         n, n, 1, 1, 1, 1, &
         fr, fi, fft_precompute_slots(slot_idx)%work_r1, fft_precompute_slots(slot_idx)%work_i1, &
         fft_precompute_slots(slot_idx)%work_r2, fft_precompute_slots(slot_idx)%work_i2, &
         sign, fft_precompute_slots(slot_idx)%stage_radix, fft_precompute_slots(slot_idx)%twc, &
         fft_precompute_slots(slot_idx)%tws)
      fr = fft_precompute_slots(slot_idx)%work_r1
      fi = fft_precompute_slots(slot_idx)%work_i1

      if (sign >= 0) then
         scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
         fr = scale * fr
         fi = scale * fi
      end if
   end subroutine ISFFT_1D_FFT_execute_reference

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

   subroutine ISFFT_1D_FFT_c2c_batch(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(:, :, :), fi(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: values(:, :, :)
      integer :: slot, entry, method, width, total, first, count, line, id, j, k

      if (size(fr, 1) /= n .or. any(shape(fr) /= shape(fi))) error stop 'split FFT shape mismatch'
      call Require_FFT_Precompute(n, slot)
      values(1:n, 1:16, 1:1) => fft_precompute_slots(slot)%split_work(1:16 * n)
      total = size(fr, 2) * size(fr, 3)
      method = 1
      width = 1
      do entry = 1, fft_precompute_slots(slot)%complex_strategy_count
         if (fft_precompute_slots(slot)%complex_strategy_lines(entry) /= total) cycle
         method = fft_precompute_slots(slot)%complex_strategy_method(entry)
         width = fft_precompute_slots(slot)%complex_strategy_width(entry)
         exit
      end do
      do first = 0, total - 1, 16
         count = min(16, total - first)
         do line = 1, count
            id = first + line - 1
            j = mod(id, size(fr, 2)) + 1
            k = id / size(fr, 2) + 1
            values(:, line, 1) = cmplx(fr(:, j, k), fi(:, j, k), ISFFT_REAL_KIND)
         end do
         call ISFFT_Execute_Complex_Plan(n, values(:, 1:count, :), sign, method, width)
         do line = 1, count
            id = first + line - 1
            j = mod(id, size(fr, 2)) + 1
            k = id / size(fr, 2) + 1
            fr(:, j, k) = real(values(:, line, 1), ISFFT_REAL_KIND)
            fi(:, j, k) = aimag(values(:, line, 1))
         end do
      end do
   end subroutine ISFFT_1D_FFT_c2c_batch

   !------------------------------------------------------------------------------
   !! Native complex butterflies and scalar cores
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_radix2_complex_stages(len, f, twiddle)
      integer, intent(in) :: len
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: f(len)
      complex(kind=ISFFT_REAL_KIND), intent(in) :: twiddle(len)
      integer :: i, ii, j, l, m, iistep, kstep
      complex(kind=ISFFT_REAL_KIND) :: value

      l = 1
      do while (l < len)
         iistep = 2 * l
         kstep = len / iistep
         do i = 1, len, iistep
            do m = 0, l - 1
               ii = i + m
               j = ii + l
               value = twiddle(1 + m * kstep) * f(j)
               f(j) = f(ii) - value
               f(ii) = f(ii) + value
            end do
         end do
         l = iistep
      end do
   end subroutine ISFFT_1D_FFT_radix2_complex_stages

   subroutine ISFFT_1D_FFT_radix4_complex_stages(len, f, twiddle, rotation_sign)
      integer, intent(in) :: len
      integer, intent(in) :: rotation_sign
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: f(len)
      complex(kind=ISFFT_REAL_KIND), intent(in) :: twiddle(len)
      integer :: i, group, offset, quarter, stage_size, twiddle_step
      complex(kind=ISFFT_REAL_KIND) :: a, b, c, d, rotated_b, rotated_d

      stage_size = 4
      do while (stage_size <= len)
         quarter = stage_size / 4
         twiddle_step = len / stage_size
         do group = 1, len, stage_size
            do offset = 0, quarter - 1
               i = group + offset
               a = f(i)
               b = twiddle(1 + offset * twiddle_step) * f(i + quarter)
               c = twiddle(1 + 2 * offset * twiddle_step) * f(i + 2 * quarter)
               d = twiddle(1 + 3 * offset * twiddle_step) * f(i + 3 * quarter)
               rotated_b = cmplx(-real(rotation_sign, ISFFT_REAL_KIND) * aimag(b), &
                  real(rotation_sign, ISFFT_REAL_KIND) * &
                  real(b, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               rotated_d = cmplx(-real(rotation_sign, ISFFT_REAL_KIND) * aimag(d), &
                  real(rotation_sign, ISFFT_REAL_KIND) * &
                  real(d, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               f(i) = a + b + c + d
               f(i + quarter) = a + rotated_b - c - rotated_d
               f(i + 2 * quarter) = a - b + c - d
               f(i + 3 * quarter) = a - rotated_b - c + rotated_d
            end do
         end do
         stage_size = 4 * stage_size
      end do
   end subroutine ISFFT_1D_FFT_radix4_complex_stages

   subroutine ISFFT_1D_FFT_core_radix4_complex(len, f, sign, slot_idx)
      integer, intent(in) :: len, sign, slot_idx
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: f(len)
      integer :: i, j
      real(kind=ISFFT_REAL_KIND) :: scale

      if (sign >= 0) then
         scale = 1.0_ISFFT_REAL_KIND / real(len, kind=ISFFT_REAL_KIND)
         do i = 1, len
            j = fft_precompute_slots(slot_idx)%rev_radix4(i)
            fft_precompute_slots(slot_idx)%work_c1(j) = scale * f(i)
         end do
      else
         do i = 1, len
            j = fft_precompute_slots(slot_idx)%rev_radix4(i)
            fft_precompute_slots(slot_idx)%work_c1(j) = f(i)
         end do
      end if
      f = fft_precompute_slots(slot_idx)%work_c1(1:len)

      if (sign >= 0) then
         call ISFFT_1D_FFT_radix4_complex_stages(len, f, fft_precompute_slots(slot_idx)%twiddle_forward, -1)
      else
         call ISFFT_1D_FFT_radix4_complex_stages(len, f, fft_precompute_slots(slot_idx)%twiddle_backward, +1)
      end if
   end subroutine ISFFT_1D_FFT_core_radix4_complex

   subroutine ISFFT_1D_FFT_radix4_batch_stages( &
      n, lines, sign, slot_idx, &
      batch_width, defer_forward_scale)
      integer, intent(in) :: n, sign, slot_idx, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: lines(:, :)
      logical, intent(in), optional :: defer_forward_scale
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      complex(kind=ISFFT_REAL_KIND) :: a, b, c, d, rotated_b, rotated_d, tw1, tw2, tw3
      complex(kind=ISFFT_REAL_KIND) :: sum_ac, diff_ac, sum_bd, diff_rot
      integer :: first_line, width, line, i, group, offset, quarter, stage_size, twiddle_step
      integer :: input_a, input_b, input_c, input_d
      real(kind=ISFFT_REAL_KIND) :: scale, sgn
      logical :: output_written, defer_scale

      work(1:16, 1:n) => fft_precompute_slots(slot_idx)%batch_work(1:(16) * (n))
      scale = 1.0_ISFFT_REAL_KIND
      sgn = real(sign, ISFFT_REAL_KIND)
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      if (sign >= 0 .and. .not. defer_scale) then
         scale = 1.0_ISFFT_REAL_KIND / real(n, kind=ISFFT_REAL_KIND)
      end if

      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         output_written = .false.
         !! Base-4 digit reversal is an involution.  Load the four inputs for
         !! the first butterfly directly from their reversed source positions,
         !! avoiding a full write/read pass over the temporary array.
         do group = 1, n, 4
            input_a = fft_precompute_slots(slot_idx)%rev_radix4(group)
            input_b = fft_precompute_slots(slot_idx)%rev_radix4(group + 1)
            input_c = fft_precompute_slots(slot_idx)%rev_radix4(group + 2)
            input_d = fft_precompute_slots(slot_idx)%rev_radix4(group + 3)
            do line = 1, width
               a = scale * lines(input_a, first_line + line - 1)
               b = scale * lines(input_b, first_line + line - 1)
               c = scale * lines(input_c, first_line + line - 1)
               d = scale * lines(input_d, first_line + line - 1)
               rotated_b = cmplx(sgn * aimag(b), -sgn * real(b, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               rotated_d = cmplx(sgn * aimag(d), -sgn * real(d, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               sum_ac = a + c
               diff_ac = a - c
               sum_bd = b + d
               diff_rot = rotated_b - rotated_d
               work(line, group) = sum_ac + sum_bd
               work(line, group + 1) = diff_ac + diff_rot
               work(line, group + 2) = sum_ac - sum_bd
               work(line, group + 3) = diff_ac - diff_rot
            end do
         end do

         stage_size = 16
         do while (stage_size <= n)
            quarter = stage_size / 4
            twiddle_step = n / stage_size
            do group = 1, n, stage_size
               do offset = 0, quarter - 1
                  i = group + offset
                  if (sign >= 0) then
                     tw1 = fft_precompute_slots(slot_idx)%twiddle_forward(1 + offset * twiddle_step)
                     tw2 = fft_precompute_slots(slot_idx)%twiddle_forward(1 + 2 * offset * twiddle_step)
                     tw3 = fft_precompute_slots(slot_idx)%twiddle_forward(1 + 3 * offset * twiddle_step)
                  else
                     tw1 = fft_precompute_slots(slot_idx)%twiddle_backward(1 + offset * twiddle_step)
                     tw2 = fft_precompute_slots(slot_idx)%twiddle_backward(1 + 2 * offset * twiddle_step)
                     tw3 = fft_precompute_slots(slot_idx)%twiddle_backward(1 + 3 * offset * twiddle_step)
                  end if
                  do line = 1, width
                     a = work(line, i)
                     if (offset == 0) then
                        b = work(line, i + quarter)
                        c = work(line, i + 2 * quarter)
                        d = work(line, i + 3 * quarter)
                     else
                        b = tw1 * work(line, i + quarter)
                        c = tw2 * work(line, i + 2 * quarter)
                        d = tw3 * work(line, i + 3 * quarter)
                     end if
                     rotated_b = cmplx(sgn * aimag(b), -sgn * real(b, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
                     rotated_d = cmplx(sgn * aimag(d), -sgn * real(d, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
                     sum_ac = a + c
                     diff_ac = a - c
                     sum_bd = b + d
                     diff_rot = rotated_b - rotated_d
                     if (stage_size == n) then
                        lines(i, first_line + line - 1) = sum_ac + sum_bd
                        lines(i + quarter, first_line + line - 1) = diff_ac + diff_rot
                        lines(i + 2 * quarter, first_line + line - 1) = sum_ac - sum_bd
                        lines(i + 3 * quarter, first_line + line - 1) = diff_ac - diff_rot
                     else
                        work(line, i) = sum_ac + sum_bd
                        work(line, i + quarter) = diff_ac + diff_rot
                        work(line, i + 2 * quarter) = sum_ac - sum_bd
                        work(line, i + 3 * quarter) = diff_ac - diff_rot
                     end if
                  end do
               end do
            end do
            if (stage_size == n) output_written = .true.
            stage_size = 4 * stage_size
         end do

         if (.not. output_written) then
            do i = 1, n
               do line = 1, width
                  lines(i, first_line + line - 1) = work(line, i)
               end do
            end do
         end if
      end do
   end subroutine ISFFT_1D_FFT_radix4_batch_stages

   subroutine ISFFT_1D_FFT_core_radix2_complex(len, f, sign, slot_idx)
      integer, intent(in) :: len, sign, slot_idx
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: f(len)
      integer :: i, j
      real(kind=ISFFT_REAL_KIND) :: scale

      if (sign >= 0) then
         scale = 1.0_ISFFT_REAL_KIND / real(len, kind=ISFFT_REAL_KIND)
         do i = 1, len
            j = fft_precompute_slots(slot_idx)%rev(i)
            fft_precompute_slots(slot_idx)%work_c1(j) = scale * f(i)
         end do
      else
         do i = 1, len
            j = fft_precompute_slots(slot_idx)%rev(i)
            fft_precompute_slots(slot_idx)%work_c1(j) = f(i)
         end do
      end if
      f = fft_precompute_slots(slot_idx)%work_c1(1:len)

      if (sign >= 0) then
         call ISFFT_1D_FFT_radix2_complex_stages(len, f, fft_precompute_slots(slot_idx)%twiddle_forward)
      else
         call ISFFT_1D_FFT_radix2_complex_stages(len, f, fft_precompute_slots(slot_idx)%twiddle_backward)
      end if
   end subroutine ISFFT_1D_FFT_core_radix2_complex

   !------------------------------------------------------------------------------
   !! Prepared C2C strategy dispatch
   !------------------------------------------------------------------------------
   subroutine ISFFT_Execute_FFT_Complex_Strategy(n, f, sign, slot_n, method, batch_width)
      integer, intent(in) :: n, sign, slot_n, method, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: lines(:, :)
      integer :: j, k
      !! Input must not alias the recursive mixed-radix scratch arrays.
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: scalar_re(:), scalar_im(:)

      scalar_re => fft_precompute_slots(slot_n)%scalar_re
      scalar_im => fft_precompute_slots(slot_n)%scalar_im
      if (method == 3 .and. fft_precompute_slots(slot_n)%is_radix4) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_radix4_batch_stages(n, lines, sign, slot_n, batch_width)
         return
      end if
      if (method == 4 .and. fft_precompute_slots(slot_n)%is_radix2) then
         call ISFFT_1D_FFT_radix2_complex_batch_contiguous(n, f, sign, slot_n, batch_width)
         return
      end if
      if (method == 5 .and. fft_precompute_slots(slot_n)%bluestein_m > 0) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_bluestein_complex_batch(n, lines, sign, slot_n, batch_width)
         return
      end if
      if (method == 8 .and. fft_precompute_slots(slot_n)%rader_ready) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_rader_complex_batch(n, lines, sign, slot_n, batch_width)
         return
      end if
      if (method == 7 .and. .not. fft_precompute_slots(slot_n)%is_radix2) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_mixed_complex_batch(n, lines, sign, slot_n, batch_width)
         return
      end if
      if (method == 9 .and. .not. fft_precompute_slots(slot_n)%is_radix2) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_mixed_complex_batch(n, lines, sign, slot_n, batch_width, .true.)
         return
      end if
      if (method == 10 .and. ISFFT_Is_Three_Times_Power_Of_Two(n)) then
         lines(1:n, 1:size(f, 2) * size(f, 3)) => f
         call ISFFT_1D_FFT_radix3x2_complex_batch(n, lines, sign, slot_n, batch_width)
         return
      end if
      do k = 1, size(f, 3)
         do j = 1, size(f, 2)
            if (method == 11 .or. method == 12) then
               scalar_re = real(f(:, j, k), ISFFT_REAL_KIND)
               scalar_im = aimag(f(:, j, k))
               if (method == 11) then
                  call ISFFT_1D_FFT_core_bluestein(n, scalar_re, scalar_im, sign, slot_n)
               else
                  call ISFFT_1D_FFT_core_rader(n, scalar_re, scalar_im, sign, slot_n)
               end if
               f(:, j, k) = cmplx(scalar_re, scalar_im, ISFFT_REAL_KIND)
            else if (method == 6 .and. .not. fft_precompute_slots(slot_n)%is_radix2) then
               scalar_re = real(f(:, j, k), kind=ISFFT_REAL_KIND)
               scalar_im = aimag(f(:, j, k))
               call ISFFT_1D_FFT_execute_reference(n, scalar_re, scalar_im, sign, slot_n)
               f(:, j, k) = cmplx(scalar_re, scalar_im, kind=ISFFT_REAL_KIND)
            else if (method == 2 .and. fft_precompute_slots(slot_n)%is_radix4) then
               call ISFFT_1D_FFT_core_radix4_complex(n, f(:, j, k), sign, slot_n)
            else if (fft_precompute_slots(slot_n)%is_radix2) then
               call ISFFT_1D_FFT_core_radix2_complex(n, f(:, j, k), sign, slot_n)
            else
               scalar_re = real(f(:, j, k), kind=ISFFT_REAL_KIND)
               scalar_im = aimag(f(:, j, k))
               call ISFFT_1D_FFT_execute(n, scalar_re, scalar_im, sign, slot_n)
               f(:, j, k) = cmplx(scalar_re, scalar_im, kind=ISFFT_REAL_KIND)
            end if
         end do
      end do
   end subroutine ISFFT_Execute_FFT_Complex_Strategy

   !------------------------------------------------------------------------------
   !! Contiguous batched C2C algorithms
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_radix2_complex_batch_contiguous( &
      n, f, sign, slot_n, &
      batch_width, defer_forward_scale)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      complex(kind=ISFFT_REAL_KIND), pointer :: lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      integer :: total_lines, first_line, width, line, j, rev_index
      real(kind=ISFFT_REAL_KIND) :: scale
      logical :: defer_scale

      work(1:16, 1:n) => fft_precompute_slots(slot_n)%batch_work(1:(16) * (n))
      total_lines = size(f, 2) * size(f, 3)
      lines(1:n, 1:total_lines) => f
      scale = 1.0_ISFFT_REAL_KIND
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      if (sign >= 0 .and. .not. defer_scale) scale = 1.0_ISFFT_REAL_KIND / real(n, kind=ISFFT_REAL_KIND)
      do first_line = 1, total_lines, batch_width
         width = min(batch_width, total_lines - first_line + 1)
         do j = 1, n
            rev_index = fft_precompute_slots(slot_n)%rev(j)
            do line = 1, width
               work(line, rev_index) = scale * lines(j, first_line + line - 1)
            end do
         end do
         if (sign >= 0) then
            call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_forward)
         else
            call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_backward)
         end if
         do j = 1, n
            do line = 1, width
               lines(j, first_line + line - 1) = work(line, j)
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_radix2_complex_batch_contiguous

   subroutine ISFFT_1D_FFT_radix2_complex_batch(n, lines, sign, slot_n, batch_width)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous :: lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      integer :: first_line, width, line, i, j
      real(kind=ISFFT_REAL_KIND) :: scale

      work(1:16, 1:n) => fft_precompute_slots(slot_n)%batch_work(1:(16) * (n))
      scale = 1.0_ISFFT_REAL_KIND
      if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, kind=ISFFT_REAL_KIND)
      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         do i = 1, n
            j = fft_precompute_slots(slot_n)%rev(i)
            do line = 1, width
               work(line, j) = scale * lines(i, first_line + line - 1)
            end do
         end do
         if (sign >= 0) then
            call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_forward)
         else
            call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_backward)
         end if
         do i = 1, n
            do line = 1, width
               lines(i, first_line + line - 1) = work(line, i)
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_radix2_complex_batch

   subroutine ISFFT_1D_FFT_radix3x2_complex_batch( &
      n, lines, sign, slot_n, &
      batch_width, defer_forward_scale)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: lines(:, :)
      logical, intent(in), optional :: defer_forward_scale
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      complex(kind=ISFFT_REAL_KIND) :: a, b, c, tw1, tw2, sum_bc, diff_bc, base, rotated_diff
      integer :: m, slot_m, first_line, width, line, t, k, rev_index, row
      real(kind=ISFFT_REAL_KIND) :: scale, root_im
      logical :: defer_scale

      work(1:48, 1:n / 3) => fft_precompute_slots(slot_n)%batch_work(1:(48) * (n / 3))
      if (.not. ISFFT_Is_Three_Times_Power_Of_Two(n)) &
         error stop 'ISFFT_1D_FFT_radix3x2_complex_batch: unsupported length'
      m = n / 3
      call Require_FFT_Precompute(m, slot_m)
      scale = 1.0_ISFFT_REAL_KIND
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      if (sign >= 0 .and. .not. defer_scale) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
      root_im = -real(sign, ISFFT_REAL_KIND) * sqrt(3.0_ISFFT_REAL_KIND) / 2.0_ISFFT_REAL_KIND

      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         do t = 0, m - 1
            rev_index = fft_precompute_slots(slot_m)%rev(t + 1)
            do line = 1, width
               do row = 0, 2
                  work(row * width + line, rev_index) = scale * lines(3 * t + row + 1, first_line + line - 1)
               end do
            end do
         end do
         if (sign >= 0) then
            call ISFFT_1D_FFT_radix2_batch_stages(m, work, 3 * width, fft_precompute_slots(slot_m)%twiddle_forward)
         else
            call ISFFT_1D_FFT_radix2_batch_stages(m, work, 3 * width, fft_precompute_slots(slot_m)%twiddle_backward)
         end if
         do k = 0, m - 1
            tw1 = cmplx(fft_precompute_slots(slot_n)%twc(k + 1), -real(sign, ISFFT_REAL_KIND) * &
               fft_precompute_slots(slot_n)%tws(k + 1), kind=ISFFT_REAL_KIND)
            tw2 = cmplx(fft_precompute_slots(slot_n)%twc(2 * k + 1), -real(sign, ISFFT_REAL_KIND) * &
               fft_precompute_slots(slot_n)%tws(2 * k + 1), kind=ISFFT_REAL_KIND)
            do line = 1, width
               a = work(line, k + 1)
               b = tw1 * work(width + line, k + 1)
               c = tw2 * work(2 * width + line, k + 1)
               sum_bc = b + c
               diff_bc = b - c
               base = a - 0.5_ISFFT_REAL_KIND * sum_bc
               rotated_diff = &
                  cmplx(-root_im * aimag(diff_bc), root_im * real(diff_bc, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               lines(k + 1, first_line + line - 1) = a + sum_bc
               lines(k + m + 1, first_line + line - 1) = base + rotated_diff
               lines(k + 2 * m + 1, first_line + line - 1) = base - rotated_diff
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_radix3x2_complex_batch

   subroutine ISFFT_1D_FFT_bluestein_complex_batch(n, lines, sign, slot_n, batch_width)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      complex(kind=ISFFT_REAL_KIND) :: swap_values(16), chirp, kernel_value
      integer :: convolution_len, first_line, width, line, i, j
      real(kind=ISFFT_REAL_KIND) :: scale

      work(1:16, 1:fft_precompute_slots(slot_n)%bluestein_m) => fft_precompute_slots(slot_n)%batch_work(1:(16) * &
         (fft_precompute_slots(slot_n)%bluestein_m))
      convolution_len = fft_precompute_slots(slot_n)%bluestein_m
      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         work(1:width, :) = cmplx(0.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
         do i = 1, n
            if (sign >= 0) then
               chirp = cmplx(fft_precompute_slots(slot_n)%chirp_c(i), -fft_precompute_slots(slot_n)%chirp_s(i), kind= &
                  ISFFT_REAL_KIND)
            else
               chirp = cmplx(fft_precompute_slots(slot_n)%chirp_c(i), fft_precompute_slots(slot_n)%chirp_s(i), kind= &
                  ISFFT_REAL_KIND)
            end if
            j = fft_precompute_slots(slot_n)%bluestein_rev(i)
            do line = 1, width
               work(line, j) = chirp * lines(i, first_line + line - 1) / real(convolution_len, ISFFT_REAL_KIND)
            end do
         end do
         call ISFFT_1D_FFT_radix2_batch_stages(convolution_len, work, width, &
            fft_precompute_slots(slot_n)%bluestein_twiddle_forward)
         do i = 1, convolution_len
            if (sign >= 0) then
               kernel_value = &
                  cmplx(fft_precompute_slots(slot_n)%kernel_pos_r(i), fft_precompute_slots(slot_n)%kernel_pos_i(i), &
                  kind=ISFFT_REAL_KIND)
            else
               kernel_value = &
                  cmplx(fft_precompute_slots(slot_n)%kernel_neg_r(i), fft_precompute_slots(slot_n)%kernel_neg_i(i), &
                  kind=ISFFT_REAL_KIND)
            end if
            work(1:width, i) = work(1:width, i) * kernel_value
         end do
         !! The frequency-domain product is in natural order.  Bit-reverse it
         !! in place for the inverse radix-2 stages instead of clearing and
         !! filling a second convolution-sized array for every line batch.
         do i = 1, convolution_len
            j = fft_precompute_slots(slot_n)%bluestein_rev(i)
            if (j > i) then
               swap_values(1:width) = work(1:width, i)
               work(1:width, i) = work(1:width, j)
               work(1:width, j) = swap_values(1:width)
            end if
         end do
         call ISFFT_1D_FFT_radix2_batch_stages(convolution_len, work, width, &
            fft_precompute_slots(slot_n)%bluestein_twiddle_backward)
         scale = 1.0_ISFFT_REAL_KIND
         if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
         do i = 1, n
            if (sign >= 0) then
               chirp = cmplx(fft_precompute_slots(slot_n)%chirp_c(i), -fft_precompute_slots(slot_n)%chirp_s(i), kind= &
                  ISFFT_REAL_KIND)
            else
               chirp = cmplx(fft_precompute_slots(slot_n)%chirp_c(i), fft_precompute_slots(slot_n)%chirp_s(i), kind= &
                  ISFFT_REAL_KIND)
            end if
            do line = 1, width
               lines(i, first_line + line - 1) = scale * chirp * work(line, i)
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_bluestein_complex_batch

   subroutine ISFFT_1D_FFT_rader_complex_batch(n, lines, sign, slot_n, batch_width)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: data(:, :)
      complex(kind=ISFFT_REAL_KIND) :: x0(16), line_sum(16), kernel_value
      integer :: q, slot_q, first_line, width, line, i, output_index
      real(kind=ISFFT_REAL_KIND) :: scale

      data(1:n - 1, 1:16) => fft_precompute_slots(slot_n)%rader_batch(1:(n - 1) * (16))
      q = n - 1
      call Require_FFT_Precompute(q, slot_q)
      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         do line = 1, width
            x0(line) = lines(1, first_line + line - 1)
            line_sum(line) = sum(lines(:, first_line + line - 1))
         end do
         do i = 1, q
            do line = 1, width
               data(i, line) = lines(fft_precompute_slots(slot_n)%rader_perm(i) + 1, first_line + line - 1)
            end do
         end do
         if (fft_precompute_slots(slot_q)%is_radix2) then
            call ISFFT_1D_FFT_radix2_complex_batch(q, data(:, 1:width), +1, slot_q, width)
         else
            call ISFFT_1D_FFT_mixed_complex_batch(q, data(:, 1:width), +1, slot_q, width)
         end if
         do i = 1, q
            if (sign >= 0) then
               kernel_value = &
                  cmplx(fft_precompute_slots(slot_n)%rader_kernel_pos_r(i), &
                  fft_precompute_slots(slot_n)%rader_kernel_pos_i(i), kind=ISFFT_REAL_KIND)
            else
               kernel_value = &
                  cmplx(fft_precompute_slots(slot_n)%rader_kernel_neg_r(i), &
                  fft_precompute_slots(slot_n)%rader_kernel_neg_i(i), kind=ISFFT_REAL_KIND)
            end if
            data(i, 1:width) = data(i, 1:width) * kernel_value
         end do
         if (fft_precompute_slots(slot_q)%is_radix2) then
            call ISFFT_1D_FFT_radix2_complex_batch(q, data(:, 1:width), -1, slot_q, width)
         else
            call ISFFT_1D_FFT_mixed_complex_batch(q, data(:, 1:width), -1, slot_q, width)
         end if
         scale = 1.0_ISFFT_REAL_KIND
         if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
         do line = 1, width
            lines(1, first_line + line - 1) = scale * line_sum(line)
         end do
         do i = 1, q
            output_index = fft_precompute_slots(slot_n)%rader_out_perm(i) + 1
            do line = 1, width
               lines(output_index, first_line + line - 1) = scale * (x0(line) + data(i, line))
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_rader_complex_batch

   subroutine ISFFT_1D_FFT_c2c_rader_batch(n, f, sign, batch_width)
      integer, intent(in) :: n, sign, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: lines(:, :)
      integer :: slot_n

      if (size(f, 1) /= n) error stop 'ISFFT_1D_FFT_c2c_rader_batch: shape mismatch'
      call Require_FFT_Precompute(n, slot_n)
      if (.not. fft_precompute_slots(slot_n)%rader_ready) &
         error stop 'ISFFT_1D_FFT_c2c_rader_batch: Rader plan unavailable'
      lines(1:n, 1:size(f, 2) * size(f, 3)) => f
      call ISFFT_1D_FFT_rader_complex_batch(n, lines, sign, slot_n, batch_width)
   end subroutine ISFFT_1D_FFT_c2c_rader_batch

   recursive subroutine ISFFT_1D_FFT_mixed_batch_recursive( &
      n_total, n_current, input_offset, input_stride, output_offset, level, &
      width, input_values, output_values, scratch_values, &
      sign, stage_radix, twc, tws)
      integer, intent(in) :: n_total, n_current, input_offset, input_stride, output_offset
      integer, intent(in) :: level, width, sign, stage_radix(:)
      complex(kind=ISFFT_REAL_KIND), intent(in) :: input_values(:, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: output_values(:, :), scratch_values(:, :)
      real(kind=ISFFT_REAL_KIND), intent(in) :: twc(:), tws(:)
      integer :: radix, child_size, child, k, j, input_index, root_index, line
      real(kind=ISFFT_REAL_KIND) :: cosine, sine
      complex(kind=ISFFT_REAL_KIND) :: even_part(16), odd_part(16), term
      integer :: index0, index1, index2, index3, index4, index5, index6
      complex(kind=ISFFT_REAL_KIND) :: twiddle, twiddle1, twiddle2, twiddle3, twiddle4, twiddle5, twiddle6

      if (n_current == 1) then
         output_values(1:width, output_offset) = input_values(1:width, input_offset)
         return
      end if
      radix = stage_radix(level)
      child_size = n_current / radix
      !! Odd prime leaves: reuse x(j)+x(r-j) and x(j)-x(r-j),
      !! then form output k and r-k together.  This halves the trigonometric
      !! products and avoids recursive calls and a full leaf copy.
      if (child_size == 1 .and. mod(radix, 2) == 1) then
         do line = 1, width
            output_values(line, output_offset) = input_values(line, input_offset)
         end do
         do j = 1, radix / 2
            do line = 1, width
               scratch_values(line, output_offset + j) = input_values(line, input_offset + j * input_stride) + &
                  input_values(line, input_offset + (radix - j) * input_stride)
               scratch_values(line, output_offset + radix - j) = input_values(line, input_offset + j * input_stride) - &
                  input_values(line, input_offset + (radix - j) * input_stride)
               output_values(line, output_offset) = output_values(line, output_offset) + &
                  scratch_values(line, output_offset + j)
            end do
         end do
         do k = 1, radix / 2
            even_part(1:width) = input_values(1:width, input_offset)
            odd_part(1:width) = cmplx(0.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
            do j = 1, radix / 2
               root_index = mod(j * k * (n_total / n_current), n_total)
               cosine = twc(root_index + 1)
               sine = real(sign, ISFFT_REAL_KIND) * tws(root_index + 1)
               do line = 1, width
                  even_part(line) = even_part(line) + cosine * scratch_values(line, output_offset + j)
                  term = scratch_values(line, output_offset + radix - j)
                  odd_part(line) = odd_part(line) + sine * &
                     cmplx(aimag(term), -real(term, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
               end do
            end do
            do line = 1, width
               output_values(line, output_offset + k) = even_part(line) + odd_part(line)
               output_values(line, output_offset + radix - k) = even_part(line) - odd_part(line)
            end do
         end do
         return
      end if
      do child = 0, radix - 1
         call ISFFT_1D_FFT_mixed_batch_recursive( &
            n_total, child_size, input_offset + child * input_stride, input_stride * radix, &
            output_offset + child * child_size, level + 1, &
            width, input_values, output_values, scratch_values, &
            sign, stage_radix, twc, tws)
      end do
      scratch_values(1:width, output_offset:output_offset + n_current - 1) = &
         output_values(1:width, output_offset:output_offset + n_current - 1)
      select case (radix)
      case (2)
         do k = 0, n_current - 1
            index0 = output_offset + mod(k, child_size)
            index1 = index0 + child_size
            root_index = mod(k * (n_total / n_current), n_total)
            twiddle1 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            output_values(1:width, output_offset + k) = scratch_values(1:width, index0) + twiddle1 * &
               scratch_values(1:width, index1)
         end do
      case (3)
         do k = 0, n_current - 1
            index0 = output_offset + mod(k, child_size)
            index1 = index0 + child_size
            index2 = index1 + child_size
            root_index = mod(k * (n_total / n_current), n_total)
            twiddle1 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(2 * k * (n_total / n_current), n_total)
            twiddle2 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            output_values(1:width, output_offset + k) = scratch_values(1:width, index0) + twiddle1 * &
               scratch_values(1:width, index1) + twiddle2 * scratch_values(1:width, index2)
         end do
      case (5)
         do k = 0, n_current - 1
            index0 = output_offset + mod(k, child_size)
            index1 = index0 + child_size
            index2 = index1 + child_size
            index3 = index2 + child_size
            index4 = index3 + child_size
            root_index = mod(k * (n_total / n_current), n_total)
            twiddle1 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(2 * k * (n_total / n_current), n_total)
            twiddle2 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(3 * k * (n_total / n_current), n_total)
            twiddle3 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(4 * k * (n_total / n_current), n_total)
            twiddle4 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            output_values(1:width, output_offset + k) = scratch_values(1:width, index0) + twiddle1 * &
               scratch_values(1:width, index1) + twiddle2 * scratch_values(1:width, index2) + twiddle3 * &
               scratch_values(1:width, index3) + twiddle4 * scratch_values(1:width, index4)
         end do
      case (7)
         do k = 0, n_current - 1
            index0 = output_offset + mod(k, child_size)
            index1 = index0 + child_size
            index2 = index1 + child_size
            index3 = index2 + child_size
            index4 = index3 + child_size
            index5 = index4 + child_size
            index6 = index5 + child_size
            root_index = mod(k * (n_total / n_current), n_total)
            twiddle1 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(2 * k * (n_total / n_current), n_total)
            twiddle2 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(3 * k * (n_total / n_current), n_total)
            twiddle3 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(4 * k * (n_total / n_current), n_total)
            twiddle4 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(5 * k * (n_total / n_current), n_total)
            twiddle5 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            root_index = mod(6 * k * (n_total / n_current), n_total)
            twiddle6 = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
               ISFFT_REAL_KIND)
            output_values(1:width, output_offset + k) = scratch_values(1:width, index0) + twiddle1 * &
               scratch_values(1:width, index1) + twiddle2 * scratch_values(1:width, index2) + twiddle3 * &
               scratch_values(1:width, index3) + twiddle4 * scratch_values(1:width, index4) + twiddle5 * &
               scratch_values(1:width, index5) + twiddle6 * scratch_values(1:width, index6)
         end do
      case default
         do k = 0, n_current - 1
            input_index = output_offset + mod(k, child_size)
            output_values(1:width, output_offset + k) = scratch_values(1:width, input_index)
            do j = 1, radix - 1
               input_index = output_offset + j * child_size + mod(k, child_size)
               root_index = mod(j * k * (n_total / n_current), n_total)
               twiddle = cmplx(twc(root_index + 1), -real(sign, ISFFT_REAL_KIND) * tws(root_index + 1), kind= &
                  ISFFT_REAL_KIND)
               output_values(1:width, output_offset + k) = output_values(1:width, output_offset + k) + twiddle * &
                  scratch_values(1:width, input_index)
            end do
         end do
      end select
   end subroutine ISFFT_1D_FFT_mixed_batch_recursive

   subroutine ISFFT_1D_FFT_mixed_complex_batch(n, lines, sign, slot_n, batch_width, reverse_order)
      integer, intent(in) :: n, sign, slot_n, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: lines(:, :)
      logical, intent(in), optional :: reverse_order
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: input_values(:, :), output_values(:, :), &
         scratch_values(:, :)
      integer :: first_line, width, line, i
      integer, pointer :: stage_order(:)
      real(kind=ISFFT_REAL_KIND) :: scale

      input_values(1:16, 1:n) => fft_precompute_slots(slot_n)%mixed_1(1:(16) * (n))
      output_values(1:16, 1:n) => fft_precompute_slots(slot_n)%mixed_2(1:(16) * (n))
      scratch_values(1:16, 1:n) => fft_precompute_slots(slot_n)%mixed_3(1:(16) * (n))
      stage_order => fft_precompute_slots(slot_n)%stage_radix
      if (present(reverse_order)) then
         if (reverse_order) stage_order => fft_precompute_slots(slot_n)%stage_radix(size(stage_order):1:-1)
      end if
      scale = 1.0_ISFFT_REAL_KIND
      if (sign >= 0) scale = 1.0_ISFFT_REAL_KIND / real(n, ISFFT_REAL_KIND)
      do first_line = 1, size(lines, 2), batch_width
         width = min(batch_width, size(lines, 2) - first_line + 1)
         do i = 1, n
            do line = 1, width
               input_values(line, i) = scale * lines(i, first_line + line - 1)
            end do
         end do
         call ISFFT_1D_FFT_mixed_batch_recursive( &
            n, n, 1, 1, 1, 1, &
            width, input_values, output_values, scratch_values, &
            sign, stage_order, fft_precompute_slots(slot_n)%twc, fft_precompute_slots(slot_n)%tws)
         do i = 1, n
            do line = 1, width
               lines(i, first_line + line - 1) = output_values(line, i)
            end do
         end do
      end do
   end subroutine ISFFT_1D_FFT_mixed_complex_batch

   !------------------------------------------------------------------------------
   !! C2C batch adapters and strided-axis execution
   !------------------------------------------------------------------------------
   !! Forward scaling may be deferred by the caller. Apply normalization once
   !! across the selected path, including its fallback.
   subroutine ISFFT_1D_FFT_c2c_complex_batch(n, f, sign, defer_forward_scale)
      integer, intent(in) :: n, sign
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      integer :: slot_n, entry, method, batch_width, total_lines
      logical :: defer_scale
      complex(kind=ISFFT_REAL_KIND), pointer :: lines(:, :)

      if (size(f, 1) /= n) error stop 'ISFFT_1D_FFT_c2c_complex_batch: shape mismatch'
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      call Require_FFT_Precompute(n, slot_n)
      total_lines = size(f, 2) * size(f, 3)
      method = 1
      batch_width = 1
      do entry = 1, fft_precompute_slots(slot_n)%complex_strategy_count
         if (fft_precompute_slots(slot_n)%complex_strategy_lines(entry) == total_lines) then
            method = fft_precompute_slots(slot_n)%complex_strategy_method(entry)
            batch_width = fft_precompute_slots(slot_n)%complex_strategy_width(entry)
            exit
         end if
      end do
      if (sign >= 0 .and. defer_scale .and. method == 3 .and. fft_precompute_slots(slot_n)%is_radix4) then
         lines(1:n, 1:total_lines) => f
         call ISFFT_1D_FFT_radix4_batch_stages(n, lines, sign, slot_n, batch_width, .true.)
      else if (sign >= 0 .and. defer_scale .and. method == 4 .and. fft_precompute_slots(slot_n)%is_radix2) then
         call ISFFT_1D_FFT_radix2_complex_batch_contiguous(n, f, sign, slot_n, batch_width, .true.)
      else if (sign >= 0 .and. defer_scale .and. method == 10 .and. ISFFT_Is_Three_Times_Power_Of_Two(n)) then
         lines(1:n, 1:total_lines) => f
         call ISFFT_1D_FFT_radix3x2_complex_batch(n, lines, sign, slot_n, batch_width, .true.)
      else
         call ISFFT_Execute_FFT_Complex_Strategy(n, f, sign, slot_n, method, batch_width)
         !! Generic, Rader and scalar candidates retain their established
         !! normalized convention.  Undo it once here so the Poisson algebra
         !! can still fold all three forward scale factors into one multiply.
         if (sign >= 0 .and. defer_scale) f = real(n, ISFFT_REAL_KIND) * f
      end if
   end subroutine ISFFT_1D_FFT_c2c_complex_batch

   subroutine ISFFT_1D_FFT_c2c_axis_batch( &
      n, f, axis, sign, &
      batch_width, defer_forward_scale, planned_method)
      integer, intent(in) :: n, axis, sign, batch_width
      complex(ISFFT_REAL_KIND), intent(inout) :: f(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      integer, intent(in), optional :: planned_method

      call ISFFT_Execute_Axis_Lines(n, f, axis, sign, batch_width, defer_forward_scale, planned_method)
   end subroutine ! ISFFT_1D_FFT_c2c_axis_batch

   !! The fixed line-index arrays and axis workspace support widths up to 16.
   !! Keep width checks, reservations and line-index storage in agreement.
   subroutine ISFFT_Execute_Axis_Lines( &
      n, f, axis, sign, &
      batch_width, defer_forward_scale, planned_method, selected_lines)
      integer, intent(in), optional :: selected_lines(:)
      integer, intent(in) :: n, axis, sign, batch_width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: f(:, :, :)
      integer, intent(in), optional :: planned_method
      logical, intent(in), optional :: defer_forward_scale
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: line_c(:)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: generic_lines(:, :)
      complex(kind=ISFFT_REAL_KIND) :: a, b, c, d, rotated_b, rotated_d, tw1, tw2, tw3
      complex(kind=ISFFT_REAL_KIND) :: sum_ac, diff_ac, sum_bd, diff_rot
      integer :: slot_n, total_lines, first_line, width, line, line_id, i, j, k, rev_index
      integer :: group, offset, quarter, stage_size, twiddle_step
      integer :: input_a, input_b, input_c, input_d
      integer :: strategy_method, strategy_width, entry
      integer :: line_i(16), line_k(16)
      real(kind=ISFFT_REAL_KIND) :: scale, sgn
      logical :: defer_scale
      logical :: output_written

      call Require_FFT_Precompute(n, slot_n)
      work(1:16, 1:n) => fft_precompute_slots(slot_n)%axis_batch(1:(16) * (n))
      line_c(1:n) => fft_precompute_slots(slot_n)%axis_lines(1:(n))
      generic_lines(1:n, 1:16) => fft_precompute_slots(slot_n)%axis_lines(1:(n) * (16))
      if (axis < 1 .or. axis > 3) error stop 'ISFFT_1D_FFT_c2c_axis_batch: invalid axis'
      if (size(f, axis) /= n) error stop 'ISFFT_1D_FFT_c2c_axis_batch: shape mismatch'
      if (batch_width < 1 .or. batch_width > 16) error stop 'ISFFT_1D_FFT_c2c_axis_batch: invalid batch width'
      !! A singleton transform is the identity; no radix-4 butterfly is valid.
      if (n == 1) return
      if (axis == 1) then
         call ISFFT_1D_FFT_c2c_complex_batch(n, f, sign)
         return
      end if

      call Require_FFT_Precompute(n, slot_n)
      total_lines = size(f, kind=ISFFT_INT64_KIND) / n
      if (present(selected_lines)) total_lines = size(selected_lines)
      strategy_method = 1
      strategy_width = batch_width
      if (.not. present(planned_method)) then
         do entry = 1, fft_precompute_slots(slot_n)%complex_strategy_count
            if (fft_precompute_slots(slot_n)%complex_strategy_lines(entry) == total_lines) then
               strategy_method = fft_precompute_slots(slot_n)%complex_strategy_method(entry)
               strategy_width = fft_precompute_slots(slot_n)%complex_strategy_width(entry)
               exit
            end if
         end do
      else
         strategy_method = planned_method
      end if
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      sgn = real(sign, ISFFT_REAL_KIND)

      !! Non-power-of-two algorithms with gather and scatter
      if (.not. fft_precompute_slots(slot_n)%is_radix2) then
         if ((strategy_method == 5 .and. fft_precompute_slots(slot_n)%bluestein_m > 0) .or. &
            (strategy_method == 8 .and. fft_precompute_slots(slot_n)%rader_ready)) then
            scale = 1.0_ISFFT_REAL_KIND
            if (sign >= 0 .and. defer_scale) scale = real(n, ISFFT_REAL_KIND)
            do first_line = 1, total_lines, strategy_width
               width = min(strategy_width, total_lines - first_line + 1)
               do j = 1, n
                  do line = 1, width
                     line_id = first_line + line - 2
                     if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
                     i = mod(line_id, size(f, 1)) + 1
                     k = line_id / size(f, 1) + 1
                     if (axis == 2) then
                        generic_lines(j, line) = f(i, j, k)
                     else
                        generic_lines(j, line) = f(i, k, j)
                     end if
                  end do
               end do
               if (strategy_method == 8) then
                  call ISFFT_1D_FFT_rader_complex_batch(n, generic_lines(:, 1:width), sign, slot_n, width)
               else
                  call ISFFT_1D_FFT_bluestein_complex_batch(n, generic_lines(:, 1:width), sign, slot_n, width)
               end if
               do j = 1, n
                  do line = 1, width
                     line_id = first_line + line - 2
                     if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
                     i = mod(line_id, size(f, 1)) + 1
                     k = line_id / size(f, 1) + 1
                     if (axis == 2) then
                        f(i, j, k) = scale * generic_lines(j, line)
                     else
                        f(i, k, j) = scale * generic_lines(j, line)
                     end if
                  end do
               end do
            end do
            return
         end if
         if (strategy_method == 7 .or. strategy_method == 9 .or. strategy_method == 10) then
            scale = 1.0_ISFFT_REAL_KIND
            if (sign >= 0 .and. defer_scale .and. strategy_method /= 10) scale = real(n, ISFFT_REAL_KIND)
            do first_line = 1, total_lines, strategy_width
               width = min(strategy_width, total_lines - first_line + 1)
               do j = 1, n
                  do line = 1, width
                     line_id = first_line + line - 2
                     if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
                     i = mod(line_id, size(f, 1)) + 1
                     k = line_id / size(f, 1) + 1
                     if (axis == 2) then
                        generic_lines(j, line) = f(i, j, k)
                     else
                        generic_lines(j, line) = f(i, k, j)
                     end if
                  end do
               end do
               if (strategy_method == 10) then
                  call ISFFT_1D_FFT_radix3x2_complex_batch(n, generic_lines(:, 1:width), sign, slot_n, width, &
                     defer_scale)
               else if (strategy_method == 9) then
                  call ISFFT_1D_FFT_mixed_complex_batch(n, generic_lines(:, 1:width), sign, slot_n, width, .true.)
               else
                  call ISFFT_1D_FFT_mixed_complex_batch(n, generic_lines(:, 1:width), sign, slot_n, width)
               end if
               do j = 1, n
                  do line = 1, width
                     line_id = first_line + line - 2
                     if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
                     i = mod(line_id, size(f, 1)) + 1
                     k = line_id / size(f, 1) + 1
                     if (axis == 2) then
                        f(i, j, k) = scale * generic_lines(j, line)
                     else
                        f(i, k, j) = scale * generic_lines(j, line)
                     end if
                  end do
               end do
            end do
            return
         end if
         if (axis == 2) then
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               k = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  fft_precompute_slots(slot_n)%scalar_re(1:n) = real(f(i, 1:n, k), kind=ISFFT_REAL_KIND)
                  fft_precompute_slots(slot_n)%scalar_im(1:n) = aimag(f(i, 1:n, k))
                  if (strategy_method == 6) then
                     call ISFFT_1D_FFT_execute_reference(n, fft_precompute_slots(slot_n)%scalar_re, &
                        fft_precompute_slots(slot_n)%scalar_im, sign, slot_n)
                  else
                     call ISFFT_1D_FFT_execute(n, fft_precompute_slots(slot_n)%scalar_re, &
                        fft_precompute_slots(slot_n)%scalar_im, sign, slot_n)
                  end if
                  scale = 1.0_ISFFT_REAL_KIND
                  if (sign >= 0 .and. defer_scale) scale = real(n, ISFFT_REAL_KIND)
                  f(i, 1:n, k) = scale * &
                     cmplx(fft_precompute_slots(slot_n)%scalar_re(1:n), fft_precompute_slots(slot_n)%scalar_im(1:n), &
                     kind=ISFFT_REAL_KIND)
               end block
            end do
         end if
         if (axis == 3) then
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               j = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  fft_precompute_slots(slot_n)%scalar_re(1:n) = real(f(i, j, 1:n), kind=ISFFT_REAL_KIND)
                  fft_precompute_slots(slot_n)%scalar_im(1:n) = aimag(f(i, j, 1:n))
                  if (strategy_method == 6) then
                     call ISFFT_1D_FFT_execute_reference(n, fft_precompute_slots(slot_n)%scalar_re, &
                        fft_precompute_slots(slot_n)%scalar_im, sign, slot_n)
                  else
                     call ISFFT_1D_FFT_execute(n, fft_precompute_slots(slot_n)%scalar_re, &
                        fft_precompute_slots(slot_n)%scalar_im, sign, slot_n)
                  end if
                  scale = 1.0_ISFFT_REAL_KIND
                  if (sign >= 0 .and. defer_scale) scale = real(n, ISFFT_REAL_KIND)
                  f(i, j, 1:n) = scale * &
                     cmplx(fft_precompute_slots(slot_n)%scalar_re(1:n), fft_precompute_slots(slot_n)%scalar_im(1:n), &
                     kind=ISFFT_REAL_KIND)
               end block
            end do
         end if
         return
      end if

      !! Power-of-two scalar and radix-2 paths
      if (.not. fft_precompute_slots(slot_n)%is_radix4 .and. batch_width == 1 .and. &
         .not. (sign >= 0 .and. defer_scale)) then
         if (axis == 2) then
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               k = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  line_c = f(i, 1:n, k)
                  call ISFFT_1D_FFT_core_radix2_complex(n, line_c, sign, slot_n)
                  f(i, 1:n, k) = line_c
               end block
            end do
         else
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               k = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  line_c = f(i, k, 1:n)
                  call ISFFT_1D_FFT_core_radix2_complex(n, line_c, sign, slot_n)
                  f(i, k, 1:n) = line_c
               end block
            end do
         end if
         return
      end if

      !! Powers of two that are not powers of four (128, 512, ...) used to
      !! fall through to one scalar transform per line.  Keep the same
      !! planner-selected width and vectorize the radix-2 butterflies across
      !! independent lines, just as the radix-4 path does.
      if (.not. fft_precompute_slots(slot_n)%is_radix4) then
         total_lines = size(f, kind=ISFFT_INT64_KIND) / n
         if (present(selected_lines)) total_lines = size(selected_lines)
         scale = 1.0_ISFFT_REAL_KIND
         if (sign >= 0 .and. .not. defer_scale) scale = 1.0_ISFFT_REAL_KIND / real(n, kind=ISFFT_REAL_KIND)
         do first_line = 1, total_lines, batch_width
            width = min(batch_width, total_lines - first_line + 1)
            do line = 1, width
               line_id = first_line + line - 2
               if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
               line_i(line) = mod(line_id, size(f, 1)) + 1
               line_k(line) = line_id / size(f, 1) + 1
            end do
            do j = 1, n
               rev_index = fft_precompute_slots(slot_n)%rev(j)
               do line = 1, width
                  if (axis == 2) then
                     work(line, rev_index) = scale * f(line_i(line), j, line_k(line))
                  else
                     work(line, rev_index) = scale * f(line_i(line), line_k(line), j)
                  end if
               end do
            end do
            if (sign >= 0) then
               call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_forward)
            else
               call ISFFT_1D_FFT_radix2_batch_stages(n, work, width, fft_precompute_slots(slot_n)%twiddle_backward)
            end if
            do j = 1, n
               do line = 1, width
                  if (axis == 2) then
                     f(line_i(line), j, line_k(line)) = work(line, j)
                  else
                     f(line_i(line), line_k(line), j) = work(line, j)
                  end if
               end do
            end do
         end do
         return
      end if

      !! Radix-4 scalar path
      if (batch_width == 1 .and. .not. (sign >= 0 .and. defer_scale)) then
         if (axis == 2) then
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               k = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  line_c = f(i, 1:n, k)
                  call ISFFT_1D_FFT_core_radix4_complex(n, line_c, sign, slot_n)
                  f(i, 1:n, k) = line_c
               end block
            end do
         else
            do line = 1, total_lines
               line_id = line - 1
               if (present(selected_lines)) line_id = selected_lines(line) - 1
               k = line_id / size(f, 1) + 1
               i = mod(line_id, size(f, 1)) + 1
               block
                  line_c = f(i, k, 1:n)
                  call ISFFT_1D_FFT_core_radix4_complex(n, line_c, sign, slot_n)
                  f(i, k, 1:n) = line_c
               end block
            end do
         end if
         return
      end if

      !! Radix-4 strided batches with a shared final scatter
      total_lines = size(f, kind=ISFFT_INT64_KIND) / n
      if (present(selected_lines)) total_lines = size(selected_lines)
      scale = 1.0_ISFFT_REAL_KIND
      if (sign >= 0 .and. .not. defer_scale) scale = 1.0_ISFFT_REAL_KIND / real(n, kind=ISFFT_REAL_KIND)
      do first_line = 1, total_lines, batch_width
         width = min(batch_width, total_lines - first_line + 1)
         output_written = .false.
         do line = 1, width
            line_id = first_line + line - 2
            if (present(selected_lines)) line_id = selected_lines(line_id + 1) - 1
            line_i(line) = mod(line_id, size(f, 1)) + 1
            line_k(line) = line_id / size(f, 1) + 1
         end do
         do group = 1, n, 4
            input_a = fft_precompute_slots(slot_n)%rev_radix4(group)
            input_b = fft_precompute_slots(slot_n)%rev_radix4(group + 1)
            input_c = fft_precompute_slots(slot_n)%rev_radix4(group + 2)
            input_d = fft_precompute_slots(slot_n)%rev_radix4(group + 3)
            do line = 1, width
               if (axis == 2) then
                  if (scale == 1.0_ISFFT_REAL_KIND) then
                     a = f(line_i(line), input_a, line_k(line))
                     b = f(line_i(line), input_b, line_k(line))
                     c = f(line_i(line), input_c, line_k(line))
                     d = f(line_i(line), input_d, line_k(line))
                  else
                     a = scale * f(line_i(line), input_a, line_k(line))
                     b = scale * f(line_i(line), input_b, line_k(line))
                     c = scale * f(line_i(line), input_c, line_k(line))
                     d = scale * f(line_i(line), input_d, line_k(line))
                  end if
               else
                  if (scale == 1.0_ISFFT_REAL_KIND) then
                     a = f(line_i(line), line_k(line), input_a)
                     b = f(line_i(line), line_k(line), input_b)
                     c = f(line_i(line), line_k(line), input_c)
                     d = f(line_i(line), line_k(line), input_d)
                  else
                     a = scale * f(line_i(line), line_k(line), input_a)
                     b = scale * f(line_i(line), line_k(line), input_b)
                     c = scale * f(line_i(line), line_k(line), input_c)
                     d = scale * f(line_i(line), line_k(line), input_d)
                  end if
               end if
               rotated_b = cmplx(sgn * aimag(b), -sgn * real(b, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               rotated_d = cmplx(sgn * aimag(d), -sgn * real(d, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
               sum_ac = a + c
               diff_ac = a - c
               sum_bd = b + d
               diff_rot = rotated_b - rotated_d
               work(line, group) = sum_ac + sum_bd
               work(line, group + 1) = diff_ac + diff_rot
               work(line, group + 2) = sum_ac - sum_bd
               work(line, group + 3) = diff_ac - diff_rot
            end do
         end do

         stage_size = 16
         do while (stage_size <= n)
            quarter = stage_size / 4
            twiddle_step = n / stage_size
            do group = 1, n, stage_size
               do offset = 0, quarter - 1
                  i = group + offset
                  if (sign >= 0) then
                     tw1 = fft_precompute_slots(slot_n)%twiddle_forward(1 + offset * twiddle_step)
                     tw2 = fft_precompute_slots(slot_n)%twiddle_forward(1 + 2 * offset * twiddle_step)
                     tw3 = fft_precompute_slots(slot_n)%twiddle_forward(1 + 3 * offset * twiddle_step)
                  else
                     tw1 = fft_precompute_slots(slot_n)%twiddle_backward(1 + offset * twiddle_step)
                     tw2 = fft_precompute_slots(slot_n)%twiddle_backward(1 + 2 * offset * twiddle_step)
                     tw3 = fft_precompute_slots(slot_n)%twiddle_backward(1 + 3 * offset * twiddle_step)
                  end if
                  do line = 1, width
                     a = work(line, i)
                     if (offset == 0) then
                        b = work(line, i + quarter)
                        c = work(line, i + 2 * quarter)
                        d = work(line, i + 3 * quarter)
                     else
                        b = tw1 * work(line, i + quarter)
                        c = tw2 * work(line, i + 2 * quarter)
                        d = tw3 * work(line, i + 3 * quarter)
                     end if
                     rotated_b = cmplx(sgn * aimag(b), -sgn * real(b, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
                     rotated_d = cmplx(sgn * aimag(d), -sgn * real(d, ISFFT_REAL_KIND), kind=ISFFT_REAL_KIND)
                     sum_ac = a + c
                     diff_ac = a - c
                     sum_bd = b + d
                     diff_rot = rotated_b - rotated_d
                     if (stage_size == n) then
                        if (axis == 2) then
                           f(line_i(line), i, line_k(line)) = sum_ac + sum_bd
                           f(line_i(line), i + quarter, line_k(line)) = diff_ac + diff_rot
                           f(line_i(line), i + 2 * quarter, line_k(line)) = sum_ac - sum_bd
                           f(line_i(line), i + 3 * quarter, line_k(line)) = diff_ac - diff_rot
                        else
                           f(line_i(line), line_k(line), i) = sum_ac + sum_bd
                           f(line_i(line), line_k(line), i + quarter) = diff_ac + diff_rot
                           f(line_i(line), line_k(line), i + 2 * quarter) = sum_ac - sum_bd
                           f(line_i(line), line_k(line), i + 3 * quarter) = diff_ac - diff_rot
                        end if
                     else
                        work(line, i) = sum_ac + sum_bd
                        work(line, i + quarter) = diff_ac + diff_rot
                        work(line, i + 2 * quarter) = sum_ac - sum_bd
                        work(line, i + 3 * quarter) = diff_ac - diff_rot
                     end if
                  end do
               end do
            end do
            if (stage_size == n) output_written = .true.
            stage_size = 4 * stage_size
         end do

         if (.not. output_written) then
            do j = 1, n
               do line = 1, width
                  if (axis == 2) then
                     f(line_i(line), j, line_k(line)) = work(line, j)
                  else
                     f(line_i(line), line_k(line), j) = work(line, j)
                  end if
               end do
            end do
         end if
      end do
   end subroutine ISFFT_Execute_Axis_Lines

   !------------------------------------------------------------------------------
   !! Real-line pairing and split-real batch butterflies
   !------------------------------------------------------------------------------
   !! Two independent real lines share one odd-length complex FFT.
   !! Match ISFFT's R2C convention: DC has unit weight, interior modes two.
   !! Caller-owned scratch avoids allocation in the execution path.
   subroutine ISFFT_1D_FFT_r2c_odd_batch(n, full_re, spectrum, packed)
      integer, intent(in) :: n
      real(kind=ISFFT_REAL_KIND), intent(in), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(out), contiguous, target :: spectrum(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous :: packed(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: input_lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: output_lines(:, :)
      complex(kind=ISFFT_REAL_KIND) :: u, v, difference
      integer :: total, pair, line, q

      total = size(full_re, 2) * size(full_re, 3)
      if (mod(n, 2) /= 1 .or. size(full_re, 1) /= n .or. size(spectrum, 1) /= n / 2 + 1 .or. size(packed, 1) /= n .or. &
         size(packed, 2) /= (total + 1) / 2 .or. size(packed, 3) /= 1 .or. size(spectrum, 2) /= size(full_re, 2) .or. &
         size(spectrum, 3) /= size(full_re, 3)) error stop 'odd real FFT: shape mismatch'
      input_lines(1:n, 1:total) => full_re
      output_lines(1:n / 2 + 1, 1:total) => spectrum
      do pair = 1, (total + 1) / 2
         line = 2 * pair - 1
         if (line < total) then
            do q = 1, n
               packed(q, pair, 1) = cmplx(input_lines(q, line), input_lines(q, line + 1), ISFFT_REAL_KIND)
            end do
         else
            packed(:, pair, 1) = cmplx(input_lines(:, line), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         end if
      end do
      call ISFFT_1D_FFT_c2c_complex_batch(n, packed, +1, .true.)
      do pair = 1, (total + 1) / 2
         line = 2 * pair - 1
         output_lines(1, line) = cmplx(real(packed(1, pair, 1), ISFFT_REAL_KIND), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         if (line < total) output_lines(1, line + 1) = &
            cmplx(aimag(packed(1, pair, 1)), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         do q = 1, n / 2
            u = packed(q + 1, pair, 1)
            v = conjg(packed(n - q + 1, pair, 1))
            output_lines(q + 1, line) = u + v
            difference = u - v
            if (line < total) output_lines(q + 1, line + 1) = &
               cmplx(aimag(difference), -real(difference, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
         end do
      end do
   end subroutine ISFFT_1D_FFT_r2c_odd_batch

   subroutine ISFFT_1D_FFT_radix2_batch_stages(n, work, width, twiddle)
      integer, intent(in) :: n, width
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: work(:, :)
      complex(kind=ISFFT_REAL_KIND), intent(in) :: twiddle(:)
      integer :: line, i, j, l, m, iistep, kstep
      complex(kind=ISFFT_REAL_KIND) :: value, tw, tw1, tw2, tw3, a, b, c, d, ac, am, bd, bm
      real(kind=ISFFT_REAL_KIND) :: rotation

      if (n >= 8) then
         !! Fixed radix-8 codelet for the first three radix-2 stages.  Besides
         !! removing loop-control overhead, it avoids every multiplication by
         !! the unity twiddle (all of stage 1, half of stage 2 and one quarter
         !! of stage 3).
         do i = 1, n, 2
            do line = 1, width
               value = work(line, i + 1)
               work(line, i + 1) = work(line, i) - value
               work(line, i) = work(line, i) + value
            end do
         end do
         tw = twiddle(1 + n / 4)
         do i = 1, n, 4
            do line = 1, width
               value = work(line, i + 2)
               work(line, i + 2) = work(line, i) - value
               work(line, i) = work(line, i) + value
               value = tw * work(line, i + 3)
               work(line, i + 3) = work(line, i + 1) - value
               work(line, i + 1) = work(line, i + 1) + value
            end do
         end do
         do i = 1, n, 8
            do line = 1, width
               value = work(line, i + 4)
               work(line, i + 4) = work(line, i) - value
               work(line, i) = work(line, i) + value
            end do
            do m = 1, 3
               j = i + m
               tw = twiddle(1 + m * n / 8)
               do line = 1, width
                  value = tw * work(line, j + 4)
                  work(line, j + 4) = work(line, j) - value
                  work(line, j) = work(line, j) + value
               end do
            end do
         end do
         l = 8
      else
         l = 1
      end if
      !! Pair adjacent binary stages. Bit-reversed input places the radix-4
      !! branches in a,c,b,d order. Keep the public normalization unchanged.
      !! This removes an intermediate work-array pass for every paired stage.
      do while (l <= n / 4)
         iistep = 4 * l
         kstep = n / iistep
         rotation = aimag(twiddle(1 + n / 4))
         do i = 1, n, iistep
            do m = 0, l - 1
               j = i + m
               tw1 = twiddle(1 + m * kstep)
               tw2 = twiddle(1 + 2 * m * kstep)
               tw3 = twiddle(1 + 3 * m * kstep)
               do line = 1, width
                  a = work(line, j)
                  if (m == 0) then
                     b = work(line, j + 2 * l)
                     c = work(line, j + l)
                     d = work(line, j + 3 * l)
                  else
                     b = tw1 * work(line, j + 2 * l)
                     c = tw2 * work(line, j + l)
                     d = tw3 * work(line, j + 3 * l)
                  end if
                  ac = a + c
                  am = a - c
                  bd = b + d
                  bm = b - d
                  value = cmplx(-rotation * aimag(bm), rotation * real(bm, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
                  work(line, j) = ac + bd
                  work(line, j + l) = am + value
                  work(line, j + 2 * l) = ac - bd
                  work(line, j + 3 * l) = am - value
               end do
            end do
         end do
         l = iistep
      end do
      do while (l < n)
         iistep = 2 * l
         kstep = n / iistep
         do i = 1, n, iistep
            j = i
            do line = 1, width
               value = work(line, j + l)
               work(line, j + l) = work(line, j) - value
               work(line, j) = work(line, j) + value
            end do
            do m = 1, l - 1
               j = i + m
               tw = twiddle(1 + m * kstep)
               do line = 1, width
                  value = tw * work(line, j + l)
                  work(line, j + l) = work(line, j) - value
                  work(line, j) = work(line, j) + value
               end do
            end do
         end do
         l = iistep
      end do
   end subroutine ISFFT_1D_FFT_radix2_batch_stages

   !------------------------------------------------------------------------------
   !! Prepared R2C execution variants
   !------------------------------------------------------------------------------
   !! Keep DC and even-length Nyquist weights separate from paired interior modes.
   !! Deferred transform scaling and deferred pair scaling are independent.
   subroutine ISFFT_Execute_R2C_Batch( &
      n, full_re, spectrum, sign, &
      batch_width, defer_forward_scale, defer_inverse_pair_scale, source)
      real(kind=ISFFT_REAL_KIND), optional, intent(in), contiguous, target :: source(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: input_re(:, :, :)
      integer, intent(in) :: n, sign, batch_width
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      logical, intent(in), optional :: defer_inverse_pair_scale
      real(kind=ISFFT_REAL_KIND), pointer :: real_lines(:, :), input_lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: spectral_lines(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: natural(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: generic(:, :)
      integer :: m, slot_n, slot_m, first_line, width, line, k, j, j_pair, total_lines
      real(kind=ISFFT_REAL_KIND) :: ur, ui, vr, vi, t1r, t1i, t2r, t2i, wr, wi, qr, qi, a0, b0, scale, recon_scale
      real(kind=ISFFT_REAL_KIND) :: inverse_pair_scale
      logical :: defer_scale, defer_pair_scale

      call Require_FFT_Precompute(n, slot_n)
      work(1:16, 1:n / 2) => fft_precompute_slots(slot_n)%r2c_pack(1:(16) * (n / 2))
      if (allocated(fft_precompute_slots(slot_n)%r2c_natural)) natural(1:16, 1:n / 2) => fft_precompute_slots(slot_n) &
         %r2c_natural(1:(16) * (n / 2))
      if (allocated(fft_precompute_slots(slot_n)%r2c_generic)) generic(1:n / 2, 1:16) => fft_precompute_slots(slot_n) &
         %r2c_generic(1:(n / 2) * (16))
      if (present(source)) then
         if (sign < 0) error stop 'R2C source is forward-only'
         input_re => source
      else
         input_re => full_re
      end if
      if (mod(n, 2) /= 0) error stop 'ISFFT_1D_FFT_r2c_complex_batch: n must be even'
      m = n / 2
      if (size(input_re, 1) /= n .or. size(spectrum, 1) /= m + 1 .or. size(input_re, 2) /= size(spectrum, 2) .or. &
         size(input_re, 3) /= size(spectrum, 3)) error stop 'ISFFT_1D_FFT_r2c_complex_batch: shape mismatch'
      call Require_FFT_Precompute(n, slot_n)
      call Require_FFT_Precompute(m, slot_m)
      input_lines(1:n, 1:size(input_re, 2) * size(input_re, 3)) => input_re
      if (sign < 0) real_lines(1:n, 1:size(input_re, 2) * size(input_re, 3)) => full_re
      spectral_lines(1:m + 1, 1:size(spectrum, 2) * size(spectrum, 3)) => spectrum
      total_lines = size(input_lines, 2)
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      defer_pair_scale = .false.
      if (present(defer_inverse_pair_scale)) defer_pair_scale = defer_inverse_pair_scale
      inverse_pair_scale = 0.5_ISFFT_REAL_KIND
      if (defer_pair_scale) inverse_pair_scale = 1.0_ISFFT_REAL_KIND

      do first_line = 1, total_lines, batch_width
         width = min(batch_width, total_lines - first_line + 1)
         if (sign >= 0) then
            recon_scale = 0.5_ISFFT_REAL_KIND
            if (fft_precompute_slots(slot_m)%is_radix2) then
               scale = 1.0_ISFFT_REAL_KIND / real(m, kind=ISFFT_REAL_KIND)
               !! In the deferred-normalization path, multiplying the packed
               !! input by two and every reconstruction term by one half is
               !! algebraically redundant.  Keep both factors at one.
               if (defer_scale .and. fft_precompute_slots(slot_m)%is_radix2) then
                  scale = 1.0_ISFFT_REAL_KIND
                  recon_scale = 1.0_ISFFT_REAL_KIND
               end if
               do k = 1, m
                  j = fft_precompute_slots(slot_m)%rev(k)
                  do line = 1, width
                     if (defer_scale) then
                        work(line, j) = &
                           cmplx(input_lines(2 * k - 1, first_line + line - 1), &
                           input_lines(2 * k, first_line + line - 1), kind= &
                           ISFFT_REAL_KIND)
                     else
                        work(line, j) = scale * &
                           cmplx(input_lines(2 * k - 1, first_line + line - 1), &
                           input_lines(2 * k, first_line + line - 1), kind= &
                           ISFFT_REAL_KIND)
                     end if
                  end do
               end do
               call ISFFT_1D_FFT_radix2_batch_stages(m, work, width, fft_precompute_slots(slot_m)%twiddle_forward)
            else
               do k = 1, m
                  do line = 1, width
                     generic(k, line) = &
                        cmplx(input_lines(2 * k - 1, first_line + line - 1), &
                        input_lines(2 * k, first_line + line - 1), kind= &
                        ISFFT_REAL_KIND)
                  end do
               end do
               if (ISFFT_Is_Three_Times_Power_Of_Two(m)) then
                  call ISFFT_1D_FFT_radix3x2_complex_batch(m, generic(:, 1:width), +1, slot_m, width)
               else
                  call ISFFT_1D_FFT_mixed_complex_batch(m, generic(:, 1:width), +1, slot_m, width)
               end if
               do k = 1, m
                  do line = 1, width
                     work(line, k) = generic(k, line)
                     if (defer_scale) work(line, k) = 2.0_ISFFT_REAL_KIND * real(m, ISFFT_REAL_KIND) * work(line, k)
                  end do
               end do
            end if
            do line = 1, width
               a0 = real(work(line, 1), kind=ISFFT_REAL_KIND)
               b0 = aimag(work(line, 1))
               if (defer_scale .and. fft_precompute_slots(slot_m)%is_radix2) then
                  spectral_lines(1, first_line + line - 1) = cmplx(a0 + b0, 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
                  spectral_lines(m + 1, first_line + line - 1) = &
                     cmplx(a0 - b0, 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
               else
                  spectral_lines(1, first_line + line - 1) = &
                     cmplx(recon_scale * (a0 + b0), 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
                  spectral_lines(m + 1, first_line + line - 1) = &
                     cmplx(recon_scale * (a0 - b0), 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
               end if
            end do
            !! k and m-k share the same two packed values.  Produce both
            !! half-spectrum entries together so the reconstruction reads and
            !! common additions are done once per conjugate pair.
            do k = 1, m / 2
               wr = fft_precompute_slots(slot_n)%twc(k + 1)
               wi = -fft_precompute_slots(slot_n)%tws(k + 1)
               do line = 1, width
                  ur = real(work(line, k + 1), kind=ISFFT_REAL_KIND)
                  ui = aimag(work(line, k + 1))
                  vr = real(work(line, m - k + 1), kind=ISFFT_REAL_KIND)
                  vi = aimag(work(line, m - k + 1))
                  if (defer_scale .and. fft_precompute_slots(slot_m)%is_radix2) then
                     t1r = ur + vr
                     t1i = ui - vi
                     t2r = ur - vr
                     t2i = ui + vi
                  else
                     t1r = recon_scale * (ur + vr)
                     t1i = recon_scale * (ui - vi)
                     t2r = recon_scale * (ur - vr)
                     t2i = recon_scale * (ui + vi)
                  end if
                  qr = wr * t2r - wi * t2i
                  qi = wr * t2i + wi * t2r
                  spectral_lines(k + 1, first_line + line - 1) = cmplx(t1r + qi, t1i - qr, kind=ISFFT_REAL_KIND)
                  spectral_lines(m - k + 1, first_line + line - 1) = cmplx(t1r - qi, -t1i - qr, kind=ISFFT_REAL_KIND)
               end do
            end do
         else
            do line = 1, width
               a0 = real(spectral_lines(1, first_line + line - 1), kind=ISFFT_REAL_KIND)
               b0 = real(spectral_lines(m + 1, first_line + line - 1), kind=ISFFT_REAL_KIND)
               if (fft_precompute_slots(slot_m)%is_radix2) then
                  work(line, fft_precompute_slots(slot_m)%rev(1)) = cmplx(a0 + b0, a0 - b0, kind=ISFFT_REAL_KIND)
               else
                  natural(line, 1) = cmplx(a0 + b0, a0 - b0, kind=ISFFT_REAL_KIND)
               end if
            end do
            !! k and m-k form one reconstruction pair.  The old loop visited
            !! both and overwrote the same two natural values twice; one
            !! representative per pair is sufficient, including the
            !! self-paired midpoint when m is even.
            do k = 1, m / 2
               wr = fft_precompute_slots(slot_n)%twc(k + 1)
               wi = fft_precompute_slots(slot_n)%tws(k + 1)
               if (fft_precompute_slots(slot_m)%is_radix2) then
                  j = fft_precompute_slots(slot_m)%rev(k + 1)
                  j_pair = fft_precompute_slots(slot_m)%rev(m - k + 1)
               end if
               do line = 1, width
                  if (defer_pair_scale) then
                     t1r = real(spectral_lines(k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND) + &
                        real(spectral_lines(m - k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND)
                     t1i = aimag(spectral_lines(k + 1, first_line + line - 1)) - &
                        aimag(spectral_lines(m - k + 1, first_line + line - 1))
                     t2r = real(spectral_lines(k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND) - &
                        real(spectral_lines(m - k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND)
                     t2i = aimag(spectral_lines(k + 1, first_line + line - 1)) + &
                        aimag(spectral_lines(m - k + 1, first_line + line - 1))
                  else
                     t1r = inverse_pair_scale * &
                        (real(spectral_lines(k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND) + &
                        real(spectral_lines(m - k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND))
                     t1i = inverse_pair_scale * &
                        (aimag(spectral_lines(k + 1, first_line + line - 1)) - &
                        aimag(spectral_lines(m - k + 1, first_line + line - 1)))
                     t2r = inverse_pair_scale * &
                        (real(spectral_lines(k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND) - &
                        real(spectral_lines(m - k + 1, first_line + line - 1), kind=ISFFT_REAL_KIND))
                     t2i = inverse_pair_scale * &
                        (aimag(spectral_lines(k + 1, first_line + line - 1)) + &
                        aimag(spectral_lines(m - k + 1, first_line + line - 1)))
                  end if
                  qr = wr * t2r - wi * t2i
                  qi = wr * t2i + wi * t2r
                  if (fft_precompute_slots(slot_m)%is_radix2) then
                     work(line, j) = cmplx(t1r - qi, t1i + qr, kind=ISFFT_REAL_KIND)
                     work(line, j_pair) = cmplx(t1r + qi, -(t1i - qr), kind=ISFFT_REAL_KIND)
                  else
                     natural(line, k + 1) = cmplx(t1r - qi, t1i + qr, kind=ISFFT_REAL_KIND)
                     natural(line, m - k + 1) = cmplx(t1r + qi, -(t1i - qr), kind=ISFFT_REAL_KIND)
                  end if
               end do
            end do
            if (fft_precompute_slots(slot_m)%is_radix2) then
               call ISFFT_1D_FFT_radix2_batch_stages(m, work, width, fft_precompute_slots(slot_m)%twiddle_backward)
            else
               do k = 1, m
                  do line = 1, width
                     generic(k, line) = natural(line, k)
                  end do
               end do
               if (ISFFT_Is_Three_Times_Power_Of_Two(m)) then
                  call ISFFT_1D_FFT_radix3x2_complex_batch(m, generic(:, 1:width), -1, slot_m, width)
               else
                  call ISFFT_1D_FFT_mixed_complex_batch(m, generic(:, 1:width), -1, slot_m, width)
               end if
               do k = 1, m
                  do line = 1, width
                     work(line, k) = generic(k, line)
                  end do
               end do
            end if
            do k = 1, m
               do line = 1, width
                  real_lines(2 * k - 1, first_line + line - 1) = real(work(line, k), kind=ISFFT_REAL_KIND)
                  real_lines(2 * k, first_line + line - 1) = aimag(work(line, k))
               end do
            end do
         end if
      end do
   end subroutine ISFFT_Execute_R2C_Batch

   subroutine ISFFT_Execute_R2C_Scalar(n, full_re, spectrum, sign, source)
      real(kind=ISFFT_REAL_KIND), optional, intent(in), contiguous, target :: source(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: input_re(:, :, :)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: spectrum(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: line_im(:)
      integer :: slot
      integer :: j, k

      if (present(source)) then
         if (sign < 0) error stop 'R2C source is forward-only'
         input_re => source
      else
         input_re => full_re
      end if
      call Require_FFT_Precompute(n, slot)
      line_im => fft_precompute_slots(slot)%r2c_im
      do k = 1, size(input_re, 3)
         do j = 1, size(input_re, 2)
            if (sign >= 0) then
               call R2C_Forward_Line(n, input_re(:, j, k), spectrum(:, j, k))
            else
               full_re(1:n / 2 + 1, j, k) = real(spectrum(:, j, k), kind=ISFFT_REAL_KIND)
               line_im(1:n / 2 + 1) = aimag(spectrum(:, j, k))
               call ISFFT_1D_FFT_r2c(n, full_re(:, j, k), line_im, -1)
            end if
         end do
      end do
   end subroutine ISFFT_Execute_R2C_Scalar

   subroutine ISFFT_Execute_R2C_Generic_Batch( &
      n, full_re, spectrum, sign, &
      batch_width, defer_forward_scale, source)
      real(kind=ISFFT_REAL_KIND), optional, intent(in), contiguous, target :: source(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: input_re(:, :, :)
      integer, intent(in) :: n, sign, batch_width
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout) :: spectrum(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: lines(:, :)
      integer :: slot_n, total_lines, first_line, width, line, line_id, i, j, k, m
      logical :: defer_scale

      call Require_FFT_Precompute(n, slot_n)
      lines(1:n, 1:16) => fft_precompute_slots(slot_n)%r2c_lines(1:(n) * (16))
      if (present(source)) then
         if (sign < 0) error stop 'R2C source is forward-only'
         input_re => source
      else
         input_re => full_re
      end if
      m = n / 2
      call Require_FFT_Precompute(n, slot_n)
      total_lines = size(input_re, 2) * size(input_re, 3)
      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      do first_line = 1, total_lines, batch_width
         width = min(batch_width, total_lines - first_line + 1)
         do line = 1, width
            line_id = first_line + line - 2
            j = mod(line_id, size(input_re, 2)) + 1
            k = line_id / size(input_re, 2) + 1
            if (sign >= 0) then
               lines(:, line) = cmplx(input_re(:, j, k), 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
            else
               lines(1:m + 1, line) = spectrum(:, j, k)
               !! Public internal R2C convention doubles paired interior modes.
               lines(2:(n + 1) / 2, line) = 0.5_ISFFT_REAL_KIND * lines(2:(n + 1) / 2, line)
               do i = m + 2, n
                  lines(i, line) = conjg(lines(n - i + 2, line))
               end do
            end if
         end do
         if (ISFFT_Is_Three_Times_Power_Of_Two(n)) then
            call ISFFT_1D_FFT_radix3x2_complex_batch(n, lines(:, 1:width), sign, slot_n, width, defer_scale)
         else
            call ISFFT_1D_FFT_bluestein_complex_batch(n, lines(:, 1:width), sign, slot_n, width)
         end if
         do line = 1, width
            line_id = first_line + line - 2
            j = mod(line_id, size(input_re, 2)) + 1
            k = line_id / size(input_re, 2) + 1
            if (sign >= 0) then
               spectrum(:, j, k) = lines(1:m + 1, line)
               spectrum(2:(n + 1) / 2, j, k) = 2.0_ISFFT_REAL_KIND * spectrum(2:(n + 1) / 2, j, k)
               if (defer_scale .and. .not. ISFFT_Is_Three_Times_Power_Of_Two(n)) spectrum(:, j, k) = &
                  real(n, ISFFT_REAL_KIND) * spectrum(:, j, k)
            else
               full_re(:, j, k) = real(lines(:, line), kind=ISFFT_REAL_KIND)
            end if
         end do
      end do
   end subroutine ISFFT_Execute_R2C_Generic_Batch

   !! Two real lines share one complex FFT. Planner-selected tiles bound scratch
   !! memory; the unpaired final line uses a zero imaginary input.
   subroutine ISFFT_Execute_R2C_Odd_Batch( &
      n, full_re, spectrum, sign, &
      defer_scale, pair_width, source)
      real(kind=ISFFT_REAL_KIND), optional, intent(in), contiguous, target :: source(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: input_re(:, :, :)
      integer, intent(in) :: n, sign
      logical, intent(in) :: defer_scale
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer :: inputs(:, :), real_output(:, :)
      complex(kind=ISFFT_REAL_KIND), pointer :: outputs(:, :)
      integer, intent(in) :: pair_width
      complex(kind=ISFFT_REAL_KIND) :: a, b, u, v
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: work(:, :, :)
      integer :: slot_n
      real(kind=ISFFT_REAL_KIND) :: scale, value_im
      integer :: total, first, width, pair, line, k, h

      call Require_FFT_Precompute(n, slot_n)
      work(1:n, 1:max(1, pair_width), 1:1) => fft_precompute_slots(slot_n)%r2c_lines(1:(n) * (max(1, pair_width)) * (1))
      if (present(source)) then
         if (sign < 0) error stop 'R2C source is forward-only'
         input_re => source
      else
         input_re => full_re
      end if
      h = n / 2 + 1
      total = size(input_re, 2) * size(input_re, 3)
      if (size(input_re, 1) /= n .or. size(spectrum, 1) /= h .or. size(spectrum, 2) /= size(input_re, 2) .or. &
         size(spectrum, 3) /= size(input_re, 3)) error stop 'ISFFT_Execute_R2C_Odd_Batch: shape mismatch'
      inputs(1:n, 1:total) => input_re
      if (sign < 0) real_output(1:n, 1:total) => full_re
      outputs(1:h, 1:total) => spectrum
      scale = 1.0_ISFFT_REAL_KIND
      if (defer_scale) scale = real(n, ISFFT_REAL_KIND)
      do first = 1, total, 2 * max(1, pair_width)
         width = min(max(1, pair_width), (total - first + 2) / 2)
         do pair = 1, width
            line = first + 2 * (pair - 1)
            if (sign >= 0) then
               do k = 1, n
                  value_im = 0.0_ISFFT_REAL_KIND
                  if (line < total) value_im = inputs(k, line + 1)
                  work(k, pair, 1) = cmplx(inputs(k, line), value_im, ISFFT_REAL_KIND)
               end do
            else
               b = cmplx(0.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
               if (line < total) b = outputs(1, line + 1)
               work(1, pair, 1) = &
                  cmplx(real(outputs(1, line), ISFFT_REAL_KIND), real(b, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
               do k = 2, h
                  a = 0.5_ISFFT_REAL_KIND * outputs(k, line)
                  b = cmplx(0.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
                  if (line < total) b = 0.5_ISFFT_REAL_KIND * outputs(k, line + 1)
                  work(k, pair, 1) = a + cmplx(-aimag(b), real(b, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
                  work(n - k + 2, pair, 1) = conjg(a) + cmplx(aimag(b), real(b, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
               end do
            end if
         end do
         call ISFFT_1D_FFT_c2c_complex_batch(n, work(:, 1:width, :), sign)
         do pair = 1, width
            line = first + 2 * (pair - 1)
            if (sign >= 0) then
               outputs(1, line) = &
                  cmplx(scale * real(work(1, pair, 1), ISFFT_REAL_KIND), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
               if (line < total) outputs(1, line + 1) = &
                  cmplx(scale * aimag(work(1, pair, 1)), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
               do k = 2, h
                  u = work(k, pair, 1)
                  v = conjg(work(n - k + 2, pair, 1))
                  outputs(k, line) = scale * (u + v)
                  b = u - v
                  if (line < total) outputs(k, line + 1) = scale * &
                     cmplx(aimag(b), -real(b, ISFFT_REAL_KIND), ISFFT_REAL_KIND)
               end do
            else
               real_output(:, line) = real(work(:, pair, 1), ISFFT_REAL_KIND)
               if (line < total) real_output(:, line + 1) = aimag(work(:, pair, 1))
            end if
         end do
      end do
   end subroutine ISFFT_Execute_R2C_Odd_Batch

   !------------------------------------------------------------------------------
   !! R2C batch strategy selection
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_r2c_complex_batch( &
      n, full_re, spectrum, sign, &
      defer_forward_scale, defer_inverse_pair_scale)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      logical, intent(in), optional :: defer_forward_scale
      logical, intent(in), optional :: defer_inverse_pair_scale
      integer :: slot_n, entry, method, batch_width, total_lines
      logical :: defer_scale, defer_pair_scale

      defer_scale = .false.
      if (present(defer_forward_scale)) defer_scale = defer_forward_scale
      defer_pair_scale = .false.
      if (present(defer_inverse_pair_scale)) defer_pair_scale = defer_inverse_pair_scale
      call Require_FFT_Precompute(n, slot_n)
      total_lines = size(full_re, 2) * size(full_re, 3)
      method = 1
      batch_width = 1
      !! Untuned direct calls retain the scalar path; Planner owns all selection.
      do entry = 1, fft_precompute_slots(slot_n)%r2c_strategy_count
         if (fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) == total_lines) then
            method = fft_precompute_slots(slot_n)%r2c_strategy_method(entry)
            batch_width = fft_precompute_slots(slot_n)%r2c_strategy_width(entry)
            exit
         end if
      end do
      if (method == 5) then
         call ISFFT_Execute_R2C_Odd_Batch(n, full_re, spectrum, sign, defer_scale, batch_width)
      else if (method == 2 .or. method == 4) then
         call ISFFT_Execute_R2C_Batch(n, full_re, spectrum, sign, batch_width, defer_scale, defer_pair_scale)
      else if (method == 3) then
         call ISFFT_Execute_R2C_Generic_Batch(n, full_re, spectrum, sign, batch_width, defer_scale)
      else
         call ISFFT_Execute_R2C_Scalar(n, full_re, spectrum, sign)
         if (sign >= 0 .and. defer_scale) spectrum = real(n, kind=ISFFT_REAL_KIND) * spectrum
      end if
   end subroutine ISFFT_1D_FFT_r2c_complex_batch

   !------------------------------------------------------------------------------
   !! Explicit C2C algorithm adapters
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c_reference(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Require_FFT_Precompute(n, slot_n)
      call ISFFT_1D_FFT_execute_reference(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_reference

   subroutine ISFFT_1D_FFT_c2c_rader(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Require_Rader_Precompute(n, slot_n)
      call ISFFT_1D_FFT_core_rader(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_rader

   subroutine ISFFT_1D_FFT_c2c_bluestein(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n

      call Require_FFT_Precompute(n, slot_n)
      call ISFFT_1D_FFT_core_bluestein(n, fr, fi, sign, slot_n)
   end subroutine ISFFT_1D_FFT_c2c_bluestein

   !------------------------------------------------------------------------------
   !! Scalar real-line forward and bidirectional adapters
   !------------------------------------------------------------------------------
   subroutine R2C_Forward_Line(n, fr, out)
      integer, intent(in) :: n
      real(kind=ISFFT_REAL_KIND), intent(in) :: fr(n)
      complex(kind=ISFFT_REAL_KIND), intent(out) :: out(n / 2 + 1)
      integer :: m, k, slot_n, slot_m
      real(kind=ISFFT_REAL_KIND) :: ur, ui, vr, vi, t1r, t1i, t2r, t2i, wr, wi, qr, qi, a0, b0

      m = n / 2
      call Require_FFT_Precompute(n, slot_n)
      if (mod(n, 2) /= 0) then
         associate (r => fft_precompute_slots(slot_n)%r2c_re, im => fft_precompute_slots(slot_n)%r2c_im)
            r = fr
            im = 0
            call ISFFT_1D_FFT_c2c(n, r, im, +1)
            out = cmplx(r(1:m + 1), im(1:m + 1), ISFFT_REAL_KIND)
            out(2:m + 1) = 2.0_ISFFT_REAL_KIND * out(2:m + 1)
            out(1) = cmplx(real(out(1), ISFFT_REAL_KIND), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         end associate
         return
      end if
      call Require_FFT_Precompute(m, slot_m)
      associate (yr => fft_precompute_slots(slot_n)%work_r1(1:m), yi => fft_precompute_slots(slot_n)%work_i1(1:m))
         do k = 1, M
            yr(k) = fr(2 * k - 1)
            yi(k) = fr(2 * k)
         end do

         call ISFFT_1D_FFT_execute(M, yr, yi, +1, slot_m)

         a0 = yr(1)
         b0 = yi(1)
         out(1) = cmplx(0.5_ISFFT_REAL_KIND * (a0 + b0), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         out(M + 1) = cmplx(0.5_ISFFT_REAL_KIND * (a0 - b0), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)

         do k = 1, M - 1
            ur = yr(k + 1)
            ui = yi(k + 1)
            vr = yr(M - k + 1)
            vi = yi(M - k + 1)

            t1r = 0.5_ISFFT_REAL_KIND * (ur + vr)
            t1i = 0.5_ISFFT_REAL_KIND * (ui - vi)
            t2r = 0.5_ISFFT_REAL_KIND * (ur - vr)
            t2i = 0.5_ISFFT_REAL_KIND * (ui + vi)

            wr = fft_precompute_slots(slot_n)%twc(k + 1)
            wi = -fft_precompute_slots(slot_n)%tws(k + 1)

            qr = wr * t2r - wi * t2i
            qi = wr * t2i + wi * t2r

            out(k + 1) = cmplx(t1r + qi, t1i - qr, ISFFT_REAL_KIND)
         end do
      end associate
   end subroutine ! R2C_Forward_Line

   subroutine ISFFT_1D_FFT_r2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)

      integer :: M, k, slot_n, slot_m
      real(kind=ISFFT_REAL_KIND) :: t1r, t1i, t2r, t2i
      real(kind=ISFFT_REAL_KIND) :: wr, wi, qr, qi
      real(kind=ISFFT_REAL_KIND) :: a0, b0

      if (sign >= 0) then
         call Require_FFT_Precompute(n, slot_n)
         call R2C_Forward_Line(n, fr, fft_precompute_slots(slot_n)%work_c1(1:n / 2 + 1))
         fr(1:n / 2 + 1) = real(fft_precompute_slots(slot_n)%work_c1(1:n / 2 + 1), ISFFT_REAL_KIND)
         fi(1:n / 2 + 1) = aimag(fft_precompute_slots(slot_n)%work_c1(1:n / 2 + 1))
         return
      end if
      if (n < 1) error stop 'FFT1D_r2c: n must be positive'
      M = n / 2
      if (mod(n, 2) /= 0) then
         !! Odd N has floor(N/2)+1 independent modes and no Nyquist bin.
         !! Keep the public full-size buffers and existing weighted convention:
         !! DC is F(0)/N, all other independent modes are 2*F(k)/N.
         fr(2:M + 1) = 0.5_ISFFT_REAL_KIND * fr(2:M + 1)
         fi(2:M + 1) = 0.5_ISFFT_REAL_KIND * fi(2:M + 1)
         fi(1) = 0.0_ISFFT_REAL_KIND
         do k = M + 2, n
            fr(k) = fr(n - k + 2)
            fi(k) = -fi(n - k + 2)
         end do
         call ISFFT_1D_FFT_c2c(n, fr, fi, -1)
         return
      end if

      call Require_FFT_Precompute(n, slot_n)
      call Require_FFT_Precompute(M, slot_m)

      associate (yr => fft_precompute_slots(slot_n)%work_r1(1:M), yi => fft_precompute_slots(slot_n)%work_i1(1:M))
         a0 = fr(1)
         b0 = fr(M + 1)
         yr(1) = a0 + b0
         yi(1) = a0 - b0

         do k = 1, M - 1
            t1r = 0.5_ISFFT_REAL_KIND * (fr(k + 1) + fr(M - k + 1))
            t1i = 0.5_ISFFT_REAL_KIND * (fi(k + 1) - fi(M - k + 1))
            t2r = 0.5_ISFFT_REAL_KIND * (fr(k + 1) - fr(M - k + 1))
            t2i = 0.5_ISFFT_REAL_KIND * (fi(k + 1) + fi(M - k + 1))

            wr = fft_precompute_slots(slot_n)%twc(k + 1)
            wi = fft_precompute_slots(slot_n)%tws(k + 1)

            qr = wr * t2r - wi * t2i
            qi = wr * t2i + wi * t2r

            t2r = -qi
            t2i = qr

            yr(k + 1) = t1r + t2r
            yi(k + 1) = t1i + t2i
            yr(M - k + 1) = t1r - t2r
            yi(M - k + 1) = -(t1i - t2i)
         end do

         call ISFFT_1D_FFT_execute(M, yr, yi, -1, slot_m)

         do k = 1, M
            fr(2 * k - 1) = yr(k)
            fr(2 * k) = yi(k)
         end do
      end associate

   end subroutine ISFFT_1D_FFT_r2c

   !------------------------------------------------------------------------------
   !! Split-real C2C adapters
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_c2c(n, fr, fi, sign)
      integer, intent(in) :: n, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(n), fi(n)
      integer :: slot_n
      complex(kind=ISFFT_REAL_KIND), pointer, contiguous :: line(:)

      call Require_FFT_Precompute(n, slot_n)
      !! DCT-I/DST-I and block transforms arrive here one line at a time.
      call Require_FFT_Precompute(n, slot_n)
      line(1:n) => fft_precompute_slots(slot_n)%line_work(1:(n))
      !! Use the same complex power-of-two kernels selected for compact batch
      !! transforms, based on their actual packed FFT length (n-1 for DCT-I,
      !! n+1 for DST-I), instead of always reusing the oldest split radix-2
      !! implementation.
      if (fft_precompute_slots(slot_n)%is_radix2) then
         line = cmplx(fr, fi, kind=ISFFT_REAL_KIND)
         if (fft_precompute_slots(slot_n)%is_radix4) then
            call ISFFT_1D_FFT_core_radix4_complex(n, line, sign, slot_n)
         else
            call ISFFT_1D_FFT_core_radix2_complex(n, line, sign, slot_n)
         end if
         fr = real(line, kind=ISFFT_REAL_KIND)
         fi = aimag(line)
      else
         call ISFFT_1D_FFT_execute(n, fr, fi, sign, slot_n)
      end if
   end subroutine ISFFT_1D_FFT_c2c

end module ISFFT_FFT_kernel
