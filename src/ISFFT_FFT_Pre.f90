module ISFFT_FFT_Pre
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use ISFFT_Buffer, only: &
      FFT_PRECOMPUTE_MAX, FFT_Precompute_Slot, &
      fft_precompute_slots, fft_precompute_count, fft_precompute_next_slot, &
      fft_precompute_use_clock
   implicit none

   private
   public :: FFT_Initialize, FFT_Finalize, Ensure_FFT_Precompute, &
      Ensure_Generic_N, Ensure_Generic_M, Ensure_Real_Transform_Precompute, Ensure_Rader_Precompute, &
      Build_Rev_Twiddle, Is_Power_Of_Two, &
      FFT_Print_Precompute_Cache

contains

   !----------------------------------------------------------------------
   logical function Is_Power_Of_Two(len)
      integer, intent(in) :: len
      integer :: value

      if (len < 1) then
         Is_Power_Of_Two = .false.
         return
      end if

      value = len
      do while (mod(value, 2) == 0)
         value = value / 2
      end do

      Is_Power_Of_Two = (value == 1)
   end function Is_Power_Of_Two
   !----------------------------------------------------------------------
   subroutine Validate_FFT_Length(len, caller)
      integer, intent(in) :: len
      character(len=*), intent(in) :: caller

      if (len < 1) then
         write(*,'(A,A)') trim(caller), ': len must be positive'
         error stop
      end if

      if (.not. Is_Power_Of_Two(len)) then
         write(*,'(A,A)') trim(caller), ': radix-2 FFT requires len to be a power of two'
         error stop
      end if
   end subroutine Validate_FFT_Length
   !----------------------------------------------------------------------
   subroutine Reset_FFT_Precompute_Slot(slot)
      type(FFT_Precompute_Slot), intent(inout) :: slot

      if (allocated(slot%rev)) deallocate(slot%rev)
      if (allocated(slot%stage_radix)) deallocate(slot%stage_radix)
      if (allocated(slot%stage_size)) deallocate(slot%stage_size)
      if (allocated(slot%stage_stride)) deallocate(slot%stage_stride)
      if (allocated(slot%twc)) deallocate(slot%twc)
      if (allocated(slot%tws)) deallocate(slot%tws)
      if (allocated(slot%real_twc)) deallocate(slot%real_twc)
      if (allocated(slot%real_tws)) deallocate(slot%real_tws)
      if (allocated(slot%work_r1)) deallocate(slot%work_r1)
      if (allocated(slot%work_i1)) deallocate(slot%work_i1)
      if (allocated(slot%work_r2)) deallocate(slot%work_r2)
      if (allocated(slot%work_i2)) deallocate(slot%work_i2)
      if (allocated(slot%aux_r)) deallocate(slot%aux_r)
      if (allocated(slot%aux_i)) deallocate(slot%aux_i)
      if (allocated(slot%bluestein_rev)) deallocate(slot%bluestein_rev)
      if (allocated(slot%bluestein_twc)) deallocate(slot%bluestein_twc)
      if (allocated(slot%bluestein_tws)) deallocate(slot%bluestein_tws)
      if (allocated(slot%chirp_c)) deallocate(slot%chirp_c)
      if (allocated(slot%chirp_s)) deallocate(slot%chirp_s)
      if (allocated(slot%kernel_pos_r)) deallocate(slot%kernel_pos_r)
      if (allocated(slot%kernel_pos_i)) deallocate(slot%kernel_pos_i)
      if (allocated(slot%kernel_neg_r)) deallocate(slot%kernel_neg_r)
      if (allocated(slot%kernel_neg_i)) deallocate(slot%kernel_neg_i)
      if (allocated(slot%bluestein_r)) deallocate(slot%bluestein_r)
      if (allocated(slot%bluestein_i)) deallocate(slot%bluestein_i)
      if (allocated(slot%bluestein_work_r)) deallocate(slot%bluestein_work_r)
      if (allocated(slot%bluestein_work_i)) deallocate(slot%bluestein_work_i)
      if (allocated(slot%rader_perm)) deallocate(slot%rader_perm)
      if (allocated(slot%rader_out_perm)) deallocate(slot%rader_out_perm)
      if (allocated(slot%rader_stage_radix)) deallocate(slot%rader_stage_radix)
      if (allocated(slot%rader_rev)) deallocate(slot%rader_rev)
      if (allocated(slot%rader_twc)) deallocate(slot%rader_twc)
      if (allocated(slot%rader_tws)) deallocate(slot%rader_tws)
      if (allocated(slot%rader_kernel_pos_r)) deallocate(slot%rader_kernel_pos_r)
      if (allocated(slot%rader_kernel_pos_i)) deallocate(slot%rader_kernel_pos_i)
      if (allocated(slot%rader_kernel_neg_r)) deallocate(slot%rader_kernel_neg_r)
      if (allocated(slot%rader_kernel_neg_i)) deallocate(slot%rader_kernel_neg_i)
      if (allocated(slot%rader_data_r)) deallocate(slot%rader_data_r)
      if (allocated(slot%rader_data_i)) deallocate(slot%rader_data_i)
      if (allocated(slot%rader_work_r1)) deallocate(slot%rader_work_r1)
      if (allocated(slot%rader_work_i1)) deallocate(slot%rader_work_i1)
      if (allocated(slot%rader_work_r2)) deallocate(slot%rader_work_r2)
      if (allocated(slot%rader_work_i2)) deallocate(slot%rader_work_i2)

      slot%len = 0
      slot%last_used = 0
      slot%active = .false.
      slot%is_radix2 = .false.
      slot%uses_fallback = .false.
      slot%use_bluestein = .false.
      slot%use_rader = .false.
      slot%rader_ready = .false.
      slot%rader_is_radix2 = .false.
      slot%num_stages = 0
      slot%bluestein_m = 0
      slot%rader_n = 0
      slot%rader_num_stages = 0
   end subroutine Reset_FFT_Precompute_Slot
   !----------------------------------------------------------------------
   subroutine Clear_FFT_Precompute_Cache()
      integer :: i

      do i = 1, FFT_PRECOMPUTE_MAX
         call Reset_FFT_Precompute_Slot(fft_precompute_slots(i))
      end do

      fft_precompute_count = 0
      fft_precompute_next_slot = 1
      fft_precompute_use_clock = 0
   end subroutine Clear_FFT_Precompute_Cache
   !----------------------------------------------------------------------
   integer function Find_FFT_Precompute_Slot(len) result(slot_idx)
      integer, intent(in) :: len
      integer :: i

      slot_idx = 0
      do i = 1, FFT_PRECOMPUTE_MAX
         if (fft_precompute_slots(i)%active .and. fft_precompute_slots(i)%len == len) then
            slot_idx = i
            return
         end if
      end do
   end function Find_FFT_Precompute_Slot
   !----------------------------------------------------------------------
   subroutine Build_Rev_Twiddle(len, rev, twc, tws)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: len
      integer, allocatable, intent(inout) :: rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: twc(:), tws(:)
      integer :: i, j, k, mm
      real(kind=ISFFT_REAL_KIND), parameter :: pi = 3.141592653589793_ISFFT_REAL_KIND

      call Validate_FFT_Length(len, 'Build_Rev_Twiddle')

      if (allocated(rev)) deallocate(rev)
      if (allocated(twc)) deallocate(twc)
      if (allocated(tws)) deallocate(tws)
      allocate(rev(len))
      allocate(twc(len/2), tws(len/2))

      j = 1
      do i = 1, len
         rev(i) = j
         mm = len/2
         do while (mm >= 1 .and. j > mm)
            j  = j - mm
            mm = mm / 2
         end do
         j = j + mm
      end do

      do k = 0, len/2 - 1
         twc(k+1) = cos((2.0_ISFFT_REAL_KIND*pi*k)/len)
         tws(k+1) = sin((2.0_ISFFT_REAL_KIND*pi*k)/len)
      end do
   end subroutine Build_Rev_Twiddle
   !----------------------------------------------------------------------
   subroutine Radix2_Transform_Unscaled(len, fr, fi, sign, rev, twc, tws, work_r, work_i)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: len, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(len), fi(len)
      integer, intent(in) :: rev(:)
      real(kind=ISFFT_REAL_KIND), intent(in) :: twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: work_r(:), work_i(:)
      integer :: i, j, l, m, iistep, kstep
      real(kind=ISFFT_REAL_KIND) :: wr, wi, tr, ti

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
         kstep = len/iistep
         do m = 0, l - 1
            wr = twc(1 + m*kstep)
            wi = -real(sign, ISFFT_REAL_KIND)*tws(1 + m*kstep)
            do i = m + 1, len, iistep
               j = i + l
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
   end subroutine Radix2_Transform_Unscaled
   !----------------------------------------------------------------------
   subroutine Build_Bluestein_Plan(len, slot)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: len
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer :: i, mirror_index, convolution_len
      real(kind=ISFFT_REAL_KIND) :: angle
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      convolution_len = 1
      do while (convolution_len < 2*len - 1)
         convolution_len = 2*convolution_len
      end do
      slot%bluestein_m = convolution_len

      call Build_Rev_Twiddle(convolution_len, slot%bluestein_rev, &
         slot%bluestein_twc, slot%bluestein_tws)
      allocate(slot%chirp_c(len), slot%chirp_s(len))
      allocate(slot%kernel_pos_r(convolution_len), slot%kernel_pos_i(convolution_len))
      allocate(slot%kernel_neg_r(convolution_len), slot%kernel_neg_i(convolution_len))
      allocate(slot%bluestein_r(convolution_len), slot%bluestein_i(convolution_len))
      allocate(slot%bluestein_work_r(convolution_len), slot%bluestein_work_i(convolution_len))

      slot%kernel_pos_r = 0.0_ISFFT_REAL_KIND
      slot%kernel_pos_i = 0.0_ISFFT_REAL_KIND
      slot%kernel_neg_r = 0.0_ISFFT_REAL_KIND
      slot%kernel_neg_i = 0.0_ISFFT_REAL_KIND
      do i = 0, len - 1
         angle = pi*real(i, ISFFT_REAL_KIND)*real(i, ISFFT_REAL_KIND) &
            / real(len, ISFFT_REAL_KIND)
         slot%chirp_c(i+1) = cos(angle)
         slot%chirp_s(i+1) = sin(angle)
         slot%kernel_pos_r(i+1) = slot%chirp_c(i+1)
         slot%kernel_pos_i(i+1) = slot%chirp_s(i+1)
         slot%kernel_neg_r(i+1) = slot%chirp_c(i+1)
         slot%kernel_neg_i(i+1) = -slot%chirp_s(i+1)
         if (i > 0) then
            mirror_index = convolution_len - i + 1
            slot%kernel_pos_r(mirror_index) = slot%chirp_c(i+1)
            slot%kernel_pos_i(mirror_index) = slot%chirp_s(i+1)
            slot%kernel_neg_r(mirror_index) = slot%chirp_c(i+1)
            slot%kernel_neg_i(mirror_index) = -slot%chirp_s(i+1)
         end if
      end do

      call Radix2_Transform_Unscaled(convolution_len, slot%kernel_pos_r, slot%kernel_pos_i, +1, &
         slot%bluestein_rev, slot%bluestein_twc, slot%bluestein_tws, &
         slot%bluestein_work_r, slot%bluestein_work_i)
      call Radix2_Transform_Unscaled(convolution_len, slot%kernel_neg_r, slot%kernel_neg_i, +1, &
         slot%bluestein_rev, slot%bluestein_twc, slot%bluestein_tws, &
         slot%bluestein_work_r, slot%bluestein_work_i)
   end subroutine Build_Bluestein_Plan
   !----------------------------------------------------------------------
   subroutine Factorize_Length(len, factors)
      integer, intent(in) :: len
      integer, allocatable, intent(out) :: factors(:)
      integer :: count, divisor, remaining, i

      count = 0
      remaining = len
      divisor = 2
      do while (remaining > 1)
         if (mod(remaining, divisor) == 0) then
            count = count + 1
            remaining = remaining / divisor
         else if (divisor == 2) then
            divisor = 3
         else
            divisor = divisor + 2
         end if
      end do

      allocate(factors(count))
      remaining = len
      divisor = 2
      do i = 1, count
         do while (mod(remaining, divisor) /= 0)
            if (divisor == 2) then
               divisor = 3
            else
               divisor = divisor + 2
            end if
         end do
         factors(i) = divisor
         remaining = remaining / divisor
      end do
   end subroutine Factorize_Length
   !----------------------------------------------------------------------
   logical function Is_Prime(len)
      integer, intent(in) :: len
      integer :: divisor

      if (len < 2) then
         Is_Prime = .false.
         return
      end if
      if (len == 2) then
         Is_Prime = .true.
         return
      end if
      if (mod(len, 2) == 0) then
         Is_Prime = .false.
         return
      end if
      divisor = 3
      do while (divisor <= len/divisor)
         if (mod(len, divisor) == 0) then
            Is_Prime = .false.
            return
         end if
         divisor = divisor + 2
      end do
      Is_Prime = .true.
   end function Is_Prime
   !----------------------------------------------------------------------
   logical function Should_Use_Rader(len)
      integer, intent(in) :: len

      Should_Use_Rader = Is_Prime(len) .and. len > 16 .and. Is_Power_Of_Two(len-1)
   end function Should_Use_Rader
   !----------------------------------------------------------------------
   integer function Modular_Power(base, exponent, modulus) result(value)
      integer, intent(in) :: base, exponent, modulus
      integer :: power
      integer(int64) :: factor64, value64, modulus64

      value64 = 1_int64
      modulus64 = int(modulus, int64)
      factor64 = modulo(int(base, int64), modulus64)
      power = exponent
      do while (power > 0)
         if (mod(power, 2) == 1) value64 = modulo(value64*factor64, modulus64)
         power = power/2
         if (power > 0) factor64 = modulo(factor64*factor64, modulus64)
      end do
      value = int(value64)
   end function Modular_Power
   !----------------------------------------------------------------------
   integer function Find_Primitive_Root(prime) result(root)
      integer, intent(in) :: prime
      integer, allocatable :: factors(:)
      integer :: candidate, i, unique_count

      call Factorize_Length(prime-1, factors)
      unique_count = 0
      do i = 1, size(factors)
         if (i == 1) then
            unique_count = unique_count + 1
            factors(unique_count) = factors(i)
         else if (factors(i) /= factors(i-1)) then
            unique_count = unique_count + 1
            factors(unique_count) = factors(i)
         end if
      end do

      root = 0
      do candidate = 2, prime - 1
         do i = 1, unique_count
            if (Modular_Power(candidate, (prime-1)/factors(i), prime) == 1) exit
         end do
         if (i > unique_count) then
            root = candidate
            exit
         end if
      end do
      deallocate(factors)
      if (root == 0) error stop 'Find_Primitive_Root: no primitive root found'
   end function Find_Primitive_Root
   !----------------------------------------------------------------------
   subroutine Build_Rader_Plan(len, slot)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: len
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer, allocatable :: factors(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: temporary_twc(:), temporary_tws(:)
      integer :: q, root, j, k, index_value
      real(kind=ISFFT_REAL_KIND) :: angle, hr, hi, wr, wi
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      if (.not. Is_Prime(len)) return
      q = len - 1
      root = Find_Primitive_Root(len)
      slot%rader_n = q
      slot%rader_is_radix2 = Is_Power_Of_Two(q)

      call Factorize_Length(q, factors)
      slot%rader_num_stages = size(factors)
      allocate(slot%rader_stage_radix(slot%rader_num_stages))
      slot%rader_stage_radix = factors
      deallocate(factors)

      allocate(slot%rader_perm(q), slot%rader_out_perm(q))
      index_value = 1
      do j = 0, q - 1
         slot%rader_perm(j+1) = index_value
         slot%rader_out_perm(q-j) = index_value
         index_value = int(modulo(int(index_value, int64)*int(root, int64), int(len, int64)))
      end do

      allocate(slot%rader_twc(q), slot%rader_tws(q))
      do k = 0, q - 1
         slot%rader_twc(k+1) = cos(2.0_ISFFT_REAL_KIND*pi*real(k, ISFFT_REAL_KIND) &
            / real(q, ISFFT_REAL_KIND))
         slot%rader_tws(k+1) = sin(2.0_ISFFT_REAL_KIND*pi*real(k, ISFFT_REAL_KIND) &
            / real(q, ISFFT_REAL_KIND))
      end do
      if (slot%rader_is_radix2) then
         call Build_Rev_Twiddle(q, slot%rader_rev, temporary_twc, temporary_tws)
         deallocate(temporary_twc, temporary_tws)
      end if

      allocate(slot%rader_kernel_pos_r(q), slot%rader_kernel_pos_i(q))
      allocate(slot%rader_kernel_neg_r(q), slot%rader_kernel_neg_i(q))
      allocate(slot%rader_data_r(q), slot%rader_data_i(q))
      allocate(slot%rader_work_r1(q), slot%rader_work_i1(q))
      allocate(slot%rader_work_r2(q), slot%rader_work_i2(q))
      slot%rader_kernel_pos_r = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_pos_i = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_neg_r = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_neg_i = 0.0_ISFFT_REAL_KIND
      if (slot%rader_is_radix2) then
         do j = 0, q - 1
            index_value = slot%rader_out_perm(j+1)
            angle = 2.0_ISFFT_REAL_KIND*pi*real(index_value, ISFFT_REAL_KIND) &
               / real(len, ISFFT_REAL_KIND)
            slot%rader_kernel_pos_r(j+1) = cos(angle)
            slot%rader_kernel_pos_i(j+1) = -sin(angle)
            slot%rader_kernel_neg_r(j+1) = cos(angle)
            slot%rader_kernel_neg_i(j+1) = sin(angle)
         end do
         call Radix2_Transform_Unscaled(q, slot%rader_kernel_pos_r, slot%rader_kernel_pos_i, +1, &
            slot%rader_rev, slot%rader_twc, slot%rader_tws, slot%rader_work_r1, slot%rader_work_i1)
         call Radix2_Transform_Unscaled(q, slot%rader_kernel_neg_r, slot%rader_kernel_neg_i, +1, &
            slot%rader_rev, slot%rader_twc, slot%rader_tws, slot%rader_work_r1, slot%rader_work_i1)
      else
         do k = 0, q - 1
            do j = 0, q - 1
               index_value = slot%rader_out_perm(j+1)
               angle = 2.0_ISFFT_REAL_KIND*pi*real(index_value, ISFFT_REAL_KIND) &
                  / real(len, ISFFT_REAL_KIND)
               hr = cos(angle)
               hi = -sin(angle)
               wr = slot%rader_twc(mod(j*k, q)+1)
               wi = -slot%rader_tws(mod(j*k, q)+1)
               slot%rader_kernel_pos_r(k+1) = slot%rader_kernel_pos_r(k+1) + hr*wr - hi*wi
               slot%rader_kernel_pos_i(k+1) = slot%rader_kernel_pos_i(k+1) + hr*wi + hi*wr
               hi = -hi
               slot%rader_kernel_neg_r(k+1) = slot%rader_kernel_neg_r(k+1) + hr*wr - hi*wi
               slot%rader_kernel_neg_i(k+1) = slot%rader_kernel_neg_i(k+1) + hr*wi + hi*wr
            end do
         end do
      end if
      slot%rader_ready = .true.
   end subroutine Build_Rader_Plan
   !----------------------------------------------------------------------
   subroutine Build_Generic_Plan(len, slot)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: len
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer, allocatable :: factors(:)
      integer :: i, current_size
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      call Factorize_Length(len, factors)
      slot%num_stages = size(factors)
      allocate(slot%stage_radix(slot%num_stages))
      allocate(slot%stage_size(slot%num_stages))
      allocate(slot%stage_stride(slot%num_stages))

      current_size = len
      do i = 1, slot%num_stages
         slot%stage_radix(i) = factors(i)
         slot%stage_size(i) = current_size
         slot%stage_stride(i) = len / current_size
         current_size = current_size / factors(i)
      end do
      slot%uses_fallback = any(slot%stage_radix > 16)
      slot%use_bluestein = slot%uses_fallback

      allocate(slot%twc(len), slot%tws(len))
      do i = 0, len - 1
         slot%twc(i+1) = cos(2.0_ISFFT_REAL_KIND*pi*real(i, ISFFT_REAL_KIND) &
            / real(len, ISFFT_REAL_KIND))
         slot%tws(i+1) = sin(2.0_ISFFT_REAL_KIND*pi*real(i, ISFFT_REAL_KIND) &
            / real(len, ISFFT_REAL_KIND))
      end do
      allocate(slot%work_r1(len), slot%work_i1(len), slot%work_r2(len), slot%work_i2(len))
      allocate(slot%aux_r(len), slot%aux_i(len))
      if (slot%uses_fallback) call Build_Bluestein_Plan(len, slot)
      if (Should_Use_Rader(len)) call Build_Rader_Plan(len, slot)
      if (slot%rader_ready) then
         slot%use_rader = .true.
         slot%use_bluestein = .false.
      end if
      deallocate(factors)
   end subroutine Build_Generic_Plan
   !----------------------------------------------------------------------
   subroutine Ensure_FFT_Precompute(len, slot_idx)
      integer, intent(in) :: len
      integer, intent(out) :: slot_idx
      integer :: i, oldest_use

      if (len < 1) error stop 'Ensure_FFT_Precompute: len must be positive'

      slot_idx = Find_FFT_Precompute_Slot(len)
      if (slot_idx > 0) then
         fft_precompute_use_clock = fft_precompute_use_clock + 1
         fft_precompute_slots(slot_idx)%last_used = fft_precompute_use_clock
         return
      end if

      if (fft_precompute_count < FFT_PRECOMPUTE_MAX) then
         slot_idx = fft_precompute_count + 1
         fft_precompute_count = slot_idx
         if (fft_precompute_count < FFT_PRECOMPUTE_MAX) then
            fft_precompute_next_slot = fft_precompute_count + 1
         else
            fft_precompute_next_slot = 1
         end if
      else
         slot_idx = 1
         oldest_use = fft_precompute_slots(1)%last_used
         do i = 2, FFT_PRECOMPUTE_MAX
            if (fft_precompute_slots(i)%last_used < oldest_use) then
               slot_idx = i
               oldest_use = fft_precompute_slots(i)%last_used
            end if
         end do
         fft_precompute_next_slot = slot_idx
      end if

      call Reset_FFT_Precompute_Slot(fft_precompute_slots(slot_idx))
      fft_precompute_slots(slot_idx)%is_radix2 = Is_Power_Of_Two(len)
      if (fft_precompute_slots(slot_idx)%is_radix2) then
         call Build_Rev_Twiddle(len, fft_precompute_slots(slot_idx)%rev, &
            fft_precompute_slots(slot_idx)%twc, fft_precompute_slots(slot_idx)%tws)
         allocate(fft_precompute_slots(slot_idx)%work_r1(len))
         allocate(fft_precompute_slots(slot_idx)%work_i1(len))
         allocate(fft_precompute_slots(slot_idx)%work_r2(len))
         allocate(fft_precompute_slots(slot_idx)%work_i2(len))
         allocate(fft_precompute_slots(slot_idx)%aux_r(len))
         allocate(fft_precompute_slots(slot_idx)%aux_i(len))
      else
         call Build_Generic_Plan(len, fft_precompute_slots(slot_idx))
      end if
      fft_precompute_slots(slot_idx)%len = len
      fft_precompute_slots(slot_idx)%active = .true.
      fft_precompute_use_clock = fft_precompute_use_clock + 1
      fft_precompute_slots(slot_idx)%last_used = fft_precompute_use_clock
      ! call FFT_Print_Precompute_Cache()
   end subroutine Ensure_FFT_Precompute
   !----------------------------------------------------------------------
   subroutine Ensure_Rader_Precompute(len, slot_idx)
      integer, intent(in) :: len
      integer, intent(out) :: slot_idx

      call Ensure_FFT_Precompute(len, slot_idx)
      if (.not. Is_Prime(len)) error stop 'Ensure_Rader_Precompute: length must be prime'
      if (.not. fft_precompute_slots(slot_idx)%rader_ready) then
         call Build_Rader_Plan(len, fft_precompute_slots(slot_idx))
      end if
   end subroutine Ensure_Rader_Precompute
   !----------------------------------------------------------------------
   subroutine Ensure_Real_Transform_Precompute(packed_len, slot_idx)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: packed_len
      integer, intent(out), optional :: slot_idx
      integer :: resolved_slot, k
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      call Ensure_FFT_Precompute(packed_len, resolved_slot)
      if (.not. allocated(fft_precompute_slots(resolved_slot)%real_twc)) then
         allocate(fft_precompute_slots(resolved_slot)%real_twc(packed_len+1))
         allocate(fft_precompute_slots(resolved_slot)%real_tws(packed_len+1))
         do k = 0, packed_len
            fft_precompute_slots(resolved_slot)%real_twc(k+1) = &
               cos(pi*real(k, ISFFT_REAL_KIND)/real(packed_len, ISFFT_REAL_KIND))
            fft_precompute_slots(resolved_slot)%real_tws(k+1) = &
               sin(pi*real(k, ISFFT_REAL_KIND)/real(packed_len, ISFFT_REAL_KIND))
         end do
      end if
      if (present(slot_idx)) slot_idx = resolved_slot
   end subroutine Ensure_Real_Transform_Precompute
   !----------------------------------------------------------------------
   subroutine Ensure_Generic_N(n, slot_idx)
      integer, intent(in) :: n
      integer, intent(out), optional :: slot_idx
      integer :: resolved_slot

      call Ensure_FFT_Precompute(n, resolved_slot)
      if (present(slot_idx)) slot_idx = resolved_slot
   end subroutine Ensure_Generic_N
   !----------------------------------------------------------------------
   subroutine Ensure_Generic_M(m, slot_idx)
      integer, intent(in) :: m
      integer, intent(out), optional :: slot_idx
      integer :: resolved_slot

      call Ensure_FFT_Precompute(m, resolved_slot)
      if (present(slot_idx)) slot_idx = resolved_slot
   end subroutine Ensure_Generic_M
   !----------------------------------------------------------------------
   subroutine Add_Unique_Length(lengths, count, len)
      integer, intent(inout) :: lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: len

      if (len < 1) return
      if (count > 0) then
         if (any(lengths(1:count) == len)) return
      end if
      if (count >= size(lengths)) error stop 'FFT_Initialize: precompute capacity is too small'
      count = count + 1
      lengths(count) = len
   end subroutine Add_Unique_Length
   !----------------------------------------------------------------------
   subroutine Add_Periodic_Lengths(lengths, count, grid_len)
      integer, intent(inout) :: lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: grid_len

      call Add_Unique_Length(lengths, count, grid_len)
      if (mod(grid_len, 2) == 0) call Add_Unique_Length(lengths, count, grid_len/2)
   end subroutine Add_Periodic_Lengths
   !----------------------------------------------------------------------
   subroutine Add_Nonperiodic_Lengths(lengths, count, grid_len)
      integer, intent(inout) :: lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: grid_len

      if (grid_len < 2) error stop 'FFT_Initialize: nonperiodic grid length must be at least two'
      call Add_Unique_Length(lengths, count, grid_len-1)
   end subroutine Add_Nonperiodic_Lengths
   !----------------------------------------------------------------------
   subroutine FFT_Initialize(nx, ny, nz, boundary_types, active_dimension)
      integer, intent(in) :: nx, ny, nz
      integer, intent(in), optional :: boundary_types(3)
      integer, intent(in), optional :: active_dimension
      integer :: grid_lengths(3), init_lengths(FFT_PRECOMPUTE_MAX)
      integer :: i, slot_idx, init_count, dimension_count

      grid_lengths = [nx, ny, nz]
      init_lengths = 0
      init_count = 0
      dimension_count = 3
      if (present(active_dimension)) dimension_count = active_dimension
      if (dimension_count < 1 .or. dimension_count > 3) then
         error stop 'FFT_Initialize: active_dimension must be 1, 2, or 3'
      end if

      call Clear_FFT_Precompute_Cache()

      do i = 1, dimension_count
         if (present(boundary_types)) then
            select case (boundary_types(i))
            case (0)
               call Add_Periodic_Lengths(init_lengths, init_count, grid_lengths(i))
            case (1, 2)
               call Add_Nonperiodic_Lengths(init_lengths, init_count, grid_lengths(i))
            case default
               error stop 'FFT_Initialize: unsupported boundary type'
            end select
         else
            call Add_Periodic_Lengths(init_lengths, init_count, grid_lengths(i))
            call Add_Nonperiodic_Lengths(init_lengths, init_count, grid_lengths(i))
         end if
      end do

      do i = 1, init_count
         call Ensure_FFT_Precompute(init_lengths(i), slot_idx)
      end do

      do i = 1, dimension_count
         if (present(boundary_types)) then
            if (boundary_types(i) == 1 .or. boundary_types(i) == 2) then
               call Ensure_Real_Transform_Precompute(grid_lengths(i)-1)
            end if
         else
            call Ensure_Real_Transform_Precompute(grid_lengths(i)-1)
         end if
      end do

   end subroutine FFT_Initialize
   !----------------------------------------------------------------------
   subroutine FFT_Finalize()
      call Clear_FFT_Precompute_Cache()
   end subroutine FFT_Finalize
   !----------------------------------------------------------------------
   subroutine FFT_Print_Precompute_Cache(unit)
      integer, intent(in), optional :: unit
      integer :: out_unit, i
      character(len=16) :: slot_state, next_marker

      out_unit = output_unit
      if (present(unit)) out_unit = unit

      write(out_unit,'(A)') '--- FFT precompute cache ---'
      write(out_unit,'(A,I0,A,I0)') 'active slots: ', fft_precompute_count, ' / ', FFT_PRECOMPUTE_MAX
      write(out_unit,'(A,I0)') 'next replacement slot: ', fft_precompute_next_slot

      do i = 1, FFT_PRECOMPUTE_MAX
         if (fft_precompute_slots(i)%active) then
            slot_state = 'active'
         else
            slot_state = 'empty'
         end if

         if (i == fft_precompute_next_slot) then
            next_marker = '<- next'
         else
            next_marker = ''
         end if

         write(out_unit,'(A,I0,A,A,A,I0,A,L1,A,L1,A,L1,1X,A)') 'slot ', i, ': ', trim(slot_state), &
            ', len=', fft_precompute_slots(i)%len, ', radix2=', fft_precompute_slots(i)%is_radix2, &
            ', fallback=', fft_precompute_slots(i)%uses_fallback, &
            ', bluestein=', fft_precompute_slots(i)%use_bluestein, trim(next_marker)
      end do
   end subroutine FFT_Print_Precompute_Cache
   !----------------------------------------------------------------------

end module ISFFT_FFT_Pre
