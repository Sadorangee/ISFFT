!------------------------------------------------------------------------------
! ISFFT_Precompute
! FFT tables, index maps and Poisson modal coefficients.
!------------------------------------------------------------------------------
module ISFFT_Precompute
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND
   use ISFFT_Discrete_Operators, only: ISFFT_Get_Scheme_Coefficients, ISFFT_Lambda_For_BC, ISFFT_Scheme_Coefficients
   use ISFFT_Plans, only: prepared_kernel_plans
   use ISFFT_Problem_State, only: Layout
   use ISFFT_Workspace, only: Destroy_Poisson_Operator, fft_precompute_count, FFT_Precompute_Slot, fft_precompute_slots, &
      fft_precompute_use_clock, Poisson_Operator, Real_Transform_Workspace, real_workspaces, &
      sample_cache
   use ISFFT_MPI_Part, only: Layout_Box
   use ISFFT_FFT_kernel, only: ISFFT_Is_Three_Times_Power_Of_Two
   use ISFFT_Allocation, only: Reserve_FFT_Axis, Reserve_FFT_Base, Reserve_R2C_Workspace
   implicit none
   private
   public :: Build_Poisson_Operator, Build_Rev_Twiddle, Ensure_FFT_Precompute, Ensure_R2C_Precompute, Ensure_Rader_Precompute
   public :: Ensure_Real_Transform_Precompute, FFT_Finalize, FFT_Initialize, Is_Power_Of_Two, Prepare_Complex_Candidates
   public :: Prepare_R2C_Candidates, Prepare_Real_Workspace

contains

   !------------------------------------------------------------------------------
   !! Numerical precomputation and prepared storage
   !------------------------------------------------------------------------------
   logical function Is_Power_Of_Two(transform_length)
      integer, intent(in) :: transform_length
      integer :: value

      if (transform_length < 1) then
         Is_Power_Of_Two = .false.
         return
      end if

      value = transform_length
      do while (mod(value, 2) == 0)
         value = value / 2
      end do

      Is_Power_Of_Two = (value == 1)
   end function Is_Power_Of_Two

   logical function Is_Power_Of_Four(transform_length)
      integer, intent(in) :: transform_length
      integer :: value

      value = transform_length
      if (value < 1) then
         Is_Power_Of_Four = .false.
         return
      end if
      do while (mod(value, 4) == 0)
         value = value / 4
      end do
      Is_Power_Of_Four = value == 1
   end function Is_Power_Of_Four

   subroutine Build_Radix4_Reverse(transform_length, reverse_index)
      integer, intent(in) :: transform_length
      integer, allocatable, intent(out) :: reverse_index(:)
      integer :: i, value, reversed, digits, digit

      digits = 0
      value = transform_length
      do while (value > 1)
         digits = digits + 1
         value = value / 4
      end do
      allocate (reverse_index(transform_length))
      do i = 0, transform_length - 1
         value = i
         reversed = 0
         do digit = 1, digits
            reversed = 4 * reversed + mod(value, 4)
            value = value / 4
         end do
         reverse_index(i + 1) = reversed + 1
      end do
   end subroutine Build_Radix4_Reverse

   subroutine Validate_FFT_Length(transform_length, caller)
      integer, intent(in) :: transform_length
      character(len=*), intent(in) :: caller

      if (transform_length < 1) then
         write (*, '(A,A)') trim(caller), ': len must be positive'
         error stop
      end if

      if (.not. Is_Power_Of_Two(transform_length)) then
         write (*, '(A,A)') trim(caller), ': radix-2 FFT requires len to be a power of two'
         error stop
      end if
   end subroutine Validate_FFT_Length

   subroutine Reset_FFT_Precompute_Slot(slot)
      type(FFT_Precompute_Slot), intent(inout) :: slot

      if (allocated(slot%rev)) deallocate (slot%rev)
      if (allocated(slot%rev_radix4)) deallocate (slot%rev_radix4)
      if (allocated(slot%stage_radix)) deallocate (slot%stage_radix)
      if (allocated(slot%stage_size)) deallocate (slot%stage_size)
      if (allocated(slot%twc)) deallocate (slot%twc)
      if (allocated(slot%tws)) deallocate (slot%tws)
      if (allocated(slot%twiddle_forward)) deallocate (slot%twiddle_forward)
      if (allocated(slot%twiddle_backward)) deallocate (slot%twiddle_backward)
      if (allocated(slot%real_twc)) deallocate (slot%real_twc)
      if (allocated(slot%real_tws)) deallocate (slot%real_tws)
      if (allocated(slot%work_r1)) deallocate (slot%work_r1)
      if (allocated(slot%work_i1)) deallocate (slot%work_i1)
      if (allocated(slot%work_r2)) deallocate (slot%work_r2)
      if (allocated(slot%work_i2)) deallocate (slot%work_i2)
      if (allocated(slot%work_c1)) deallocate (slot%work_c1)
      if (allocated(slot%split_work)) deallocate (slot%split_work)
      if (allocated(slot%line_work)) deallocate (slot%line_work)
      if (allocated(slot%batch_work)) deallocate (slot%batch_work)
      if (allocated(slot%mixed_1)) deallocate (slot%mixed_1)
      if (allocated(slot%mixed_2)) deallocate (slot%mixed_2)
      if (allocated(slot%mixed_3)) deallocate (slot%mixed_3)
      if (allocated(slot%rader_batch)) deallocate (slot%rader_batch)
      if (allocated(slot%axis_batch)) deallocate (slot%axis_batch)
      if (allocated(slot%axis_lines)) deallocate (slot%axis_lines)
      if (allocated(slot%r2c_pack)) deallocate (slot%r2c_pack)
      if (allocated(slot%r2c_natural)) deallocate (slot%r2c_natural)
      if (allocated(slot%r2c_generic)) deallocate (slot%r2c_generic)
      if (allocated(slot%r2c_lines)) deallocate (slot%r2c_lines)
      if (allocated(slot%scalar_re)) deallocate (slot%scalar_re)
      if (allocated(slot%scalar_im)) deallocate (slot%scalar_im)
      if (allocated(slot%r2c_re)) deallocate (slot%r2c_re, slot%r2c_im)
      if (allocated(slot%aux_r)) deallocate (slot%aux_r)
      if (allocated(slot%aux_i)) deallocate (slot%aux_i)
      if (allocated(slot%bluestein_rev)) deallocate (slot%bluestein_rev)
      if (allocated(slot%bluestein_twc)) deallocate (slot%bluestein_twc)
      if (allocated(slot%bluestein_tws)) deallocate (slot%bluestein_tws)
      if (allocated(slot%bluestein_twiddle_forward)) deallocate (slot%bluestein_twiddle_forward)
      if (allocated(slot%bluestein_twiddle_backward)) deallocate (slot%bluestein_twiddle_backward)
      if (allocated(slot%chirp_c)) deallocate (slot%chirp_c)
      if (allocated(slot%chirp_s)) deallocate (slot%chirp_s)
      if (allocated(slot%kernel_pos_r)) deallocate (slot%kernel_pos_r)
      if (allocated(slot%kernel_pos_i)) deallocate (slot%kernel_pos_i)
      if (allocated(slot%kernel_neg_r)) deallocate (slot%kernel_neg_r)
      if (allocated(slot%kernel_neg_i)) deallocate (slot%kernel_neg_i)
      if (allocated(slot%bluestein_r)) deallocate (slot%bluestein_r)
      if (allocated(slot%bluestein_i)) deallocate (slot%bluestein_i)
      if (allocated(slot%bluestein_work_r)) deallocate (slot%bluestein_work_r)
      if (allocated(slot%bluestein_work_i)) deallocate (slot%bluestein_work_i)
      if (allocated(slot%rader_perm)) deallocate (slot%rader_perm)
      if (allocated(slot%rader_out_perm)) deallocate (slot%rader_out_perm)
      if (allocated(slot%rader_stage_radix)) deallocate (slot%rader_stage_radix)
      if (allocated(slot%rader_rev)) deallocate (slot%rader_rev)
      if (allocated(slot%rader_twc)) deallocate (slot%rader_twc)
      if (allocated(slot%rader_tws)) deallocate (slot%rader_tws)
      if (allocated(slot%rader_kernel_pos_r)) deallocate (slot%rader_kernel_pos_r)
      if (allocated(slot%rader_kernel_pos_i)) deallocate (slot%rader_kernel_pos_i)
      if (allocated(slot%rader_kernel_neg_r)) deallocate (slot%rader_kernel_neg_r)
      if (allocated(slot%rader_kernel_neg_i)) deallocate (slot%rader_kernel_neg_i)
      if (allocated(slot%rader_data_r)) deallocate (slot%rader_data_r)
      if (allocated(slot%rader_data_i)) deallocate (slot%rader_data_i)
      if (allocated(slot%rader_work_r1)) deallocate (slot%rader_work_r1)
      if (allocated(slot%rader_work_i1)) deallocate (slot%rader_work_i1)
      if (allocated(slot%rader_work_r2)) deallocate (slot%rader_work_r2)
      if (allocated(slot%rader_work_i2)) deallocate (slot%rader_work_i2)

      slot%transform_length = 0
      slot%last_used = 0
      slot%active = .false.
      slot%is_radix2 = .false.
      slot%complex_ready = .false.
      slot%batch_ready = .false.
      slot%candidates_ready = .false.
      slot%is_radix4 = .false.
      slot%uses_fallback = .false.
      slot%use_bluestein = .false.
      slot%use_rader = .false.
      slot%rader_ready = .false.
      slot%rader_is_radix2 = .false.
      slot%num_stages = 0
      slot%bluestein_m = 0
      slot%rader_n = 0
      slot%rader_num_stages = 0
      slot%complex_strategy_count = 0
      if (allocated(slot%complex_strategy_lines)) deallocate (slot%complex_strategy_lines)
      if (allocated(slot%complex_strategy_method)) deallocate (slot%complex_strategy_method)
      if (allocated(slot%complex_strategy_width)) deallocate (slot%complex_strategy_width)
      slot%r2c_strategy_count = 0
      if (allocated(slot%r2c_strategy_lines)) deallocate (slot%r2c_strategy_lines)
      if (allocated(slot%r2c_strategy_method)) deallocate (slot%r2c_strategy_method)
      if (allocated(slot%r2c_strategy_width)) deallocate (slot%r2c_strategy_width)
   end subroutine Reset_FFT_Precompute_Slot

   subroutine Clear_FFT_Precompute_Cache()
      integer :: i

      if (allocated(fft_precompute_slots)) then
         do i = 1, size(fft_precompute_slots)
            call Reset_FFT_Precompute_Slot(fft_precompute_slots(i))
         end do
         deallocate (fft_precompute_slots)
      end if

      if (allocated(sample_cache)) deallocate (sample_cache)
      fft_precompute_count = 0
      fft_precompute_use_clock = 0
   end subroutine Clear_FFT_Precompute_Cache

   integer function Find_FFT_Precompute_Slot(transform_length) result(slot_idx)
      integer, intent(in) :: transform_length
      integer :: i

      slot_idx = 0
      do i = 1, fft_precompute_count
         if (fft_precompute_slots(i)%active .and. fft_precompute_slots(i)%transform_length == transform_length) then
            slot_idx = i
            return
         end if
      end do
   end function Find_FFT_Precompute_Slot

   subroutine Build_Rev_Twiddle(transform_length, rev, twc, tws)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: transform_length
      integer, allocatable, intent(inout) :: rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: twc(:), tws(:)
      integer :: k
      real(kind=ISFFT_REAL_KIND), parameter :: pi = 3.141592653589793_ISFFT_REAL_KIND

      call Validate_FFT_Length(transform_length, 'Build_Rev_Twiddle')

      if (allocated(twc)) deallocate (twc)
      if (allocated(tws)) deallocate (tws)
      allocate (twc(transform_length / 2), tws(transform_length / 2))

      call Build_Bit_Reverse(transform_length, rev)

      do k = 0, transform_length / 2 - 1
         twc(k + 1) = cos((2.0_ISFFT_REAL_KIND * pi * k) / transform_length)
         tws(k + 1) = sin((2.0_ISFFT_REAL_KIND * pi * k) / transform_length)
      end do
   end subroutine Build_Rev_Twiddle

   subroutine Build_Bit_Reverse(transform_length, rev)
      integer, intent(in) :: transform_length
      integer, allocatable, intent(inout) :: rev(:)
      integer :: i, j, mm

      if (allocated(rev)) then
         if (size(rev) == transform_length) return
         deallocate (rev)
      end if
      allocate (rev(transform_length))
      j = 1
      do i = 1, transform_length
         rev(i) = j
         mm = transform_length / 2
         do while (mm >= 1 .and. j > mm)
            j = j - mm
            mm = mm / 2
         end do
         j = j + mm
      end do

   end subroutine Build_Bit_Reverse

   subroutine Radix2_Transform_Unscaled(transform_length, fr, fi, sign, rev, twc, tws, work_r, work_i)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: transform_length, sign
      real(kind=ISFFT_REAL_KIND), intent(inout) :: fr(transform_length), fi(transform_length)
      integer, intent(in) :: rev(:)
      real(kind=ISFFT_REAL_KIND), intent(in) :: twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: work_r(:), work_i(:)
      integer :: i, j, l, m, iistep, kstep
      real(kind=ISFFT_REAL_KIND) :: wr, wi, tr, ti

      do i = 1, transform_length
         j = rev(i)
         work_r(j) = fr(i)
         work_i(j) = fi(i)
      end do
      fr = work_r(1:transform_length)
      fi = work_i(1:transform_length)

      l = 1
      do while (l < transform_length)
         iistep = 2 * l
         kstep = transform_length / iistep
         do m = 0, l - 1
            wr = twc(1 + m * kstep)
            wi = -real(sign, ISFFT_REAL_KIND) * tws(1 + m * kstep)
            do i = m + 1, transform_length, iistep
               j = i + l
               tr = wr * fr(j) - wi * fi(j)
               ti = wr * fi(j) + wi * fr(j)
               fr(j) = fr(i) - tr
               fi(j) = fi(i) - ti
               fr(i) = fr(i) + tr
               fi(i) = fi(i) + ti
            end do
         end do
         l = iistep
      end do
   end subroutine Radix2_Transform_Unscaled

   subroutine Build_Bluestein_Plan(transform_length, slot)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: transform_length
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer :: i, mirror_index, convolution_len
      real(kind=ISFFT_REAL_KIND) :: angle
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      convolution_len = 1
      do while (convolution_len < 2 * transform_length - 1)
         convolution_len = 2 * convolution_len
      end do
      slot%bluestein_m = convolution_len

      call Build_Rev_Twiddle(convolution_len, slot%bluestein_rev, slot%bluestein_twc, slot%bluestein_tws)
      allocate (slot%bluestein_twiddle_forward(convolution_len))
      allocate (slot%bluestein_twiddle_backward(convolution_len))
      slot%bluestein_twiddle_forward(1:convolution_len / 2) = &
         cmplx(slot%bluestein_twc, -slot%bluestein_tws, kind=ISFFT_REAL_KIND)
      slot%bluestein_twiddle_backward(1:convolution_len / 2) = &
         cmplx(slot%bluestein_twc, slot%bluestein_tws, kind=ISFFT_REAL_KIND)
      slot%bluestein_twiddle_forward(convolution_len / 2 + 1:convolution_len) = - &
         slot%bluestein_twiddle_forward(1:convolution_len / 2)
      slot%bluestein_twiddle_backward(convolution_len / 2 + 1:convolution_len) = - &
         slot%bluestein_twiddle_backward(1:convolution_len / 2)
      allocate (slot%chirp_c(transform_length), slot%chirp_s(transform_length))
      allocate (slot%kernel_pos_r(convolution_len), slot%kernel_pos_i(convolution_len))
      allocate (slot%kernel_neg_r(convolution_len), slot%kernel_neg_i(convolution_len))
      allocate (slot%bluestein_r(convolution_len), slot%bluestein_i(convolution_len))
      allocate (slot%bluestein_work_r(convolution_len), slot%bluestein_work_i(convolution_len))

      slot%kernel_pos_r = 0.0_ISFFT_REAL_KIND
      slot%kernel_pos_i = 0.0_ISFFT_REAL_KIND
      slot%kernel_neg_r = 0.0_ISFFT_REAL_KIND
      slot%kernel_neg_i = 0.0_ISFFT_REAL_KIND
      do i = 0, transform_length - 1
         angle = pi * real(i, ISFFT_REAL_KIND) * real(i, ISFFT_REAL_KIND) / real(transform_length, ISFFT_REAL_KIND)
         slot%chirp_c(i + 1) = cos(angle)
         slot%chirp_s(i + 1) = sin(angle)
         slot%kernel_pos_r(i + 1) = slot%chirp_c(i + 1)
         slot%kernel_pos_i(i + 1) = slot%chirp_s(i + 1)
         slot%kernel_neg_r(i + 1) = slot%chirp_c(i + 1)
         slot%kernel_neg_i(i + 1) = -slot%chirp_s(i + 1)
         if (i > 0) then
            mirror_index = convolution_len - i + 1
            slot%kernel_pos_r(mirror_index) = slot%chirp_c(i + 1)
            slot%kernel_pos_i(mirror_index) = slot%chirp_s(i + 1)
            slot%kernel_neg_r(mirror_index) = slot%chirp_c(i + 1)
            slot%kernel_neg_i(mirror_index) = -slot%chirp_s(i + 1)
         end if
      end do

      call Radix2_Transform_Unscaled(convolution_len, slot%kernel_pos_r, slot%kernel_pos_i, +1, slot%bluestein_rev, &
         slot%bluestein_twc, slot%bluestein_tws, slot%bluestein_work_r, slot%bluestein_work_i)
      call Radix2_Transform_Unscaled(convolution_len, slot%kernel_neg_r, slot%kernel_neg_i, +1, slot%bluestein_rev, &
         slot%bluestein_twc, slot%bluestein_tws, slot%bluestein_work_r, slot%bluestein_work_i)
   end subroutine Build_Bluestein_Plan

   subroutine Factorize_Length(transform_length, factors)
      integer, intent(in) :: transform_length
      integer, allocatable, intent(out) :: factors(:)
      integer :: count, divisor, remaining, i

      count = 0
      remaining = transform_length
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

      allocate (factors(count))
      remaining = transform_length
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

   logical function Is_Prime(transform_length)
      integer, intent(in) :: transform_length
      integer :: divisor

      if (transform_length < 2) then
         Is_Prime = .false.
         return
      end if
      if (transform_length == 2) then
         Is_Prime = .true.
         return
      end if
      if (mod(transform_length, 2) == 0) then
         Is_Prime = .false.
         return
      end if
      divisor = 3
      do while (divisor <= transform_length / divisor)
         if (mod(transform_length, divisor) == 0) then
            Is_Prime = .false.
            return
         end if
         divisor = divisor + 2
      end do
      Is_Prime = .true.
   end function Is_Prime

   logical function Should_Use_Rader(transform_length)
      integer, intent(in) :: transform_length
      integer :: remaining

      if (.not. Is_Prime(transform_length)) then
         Should_Use_Rader = .false.
         return
      end if
      remaining = transform_length - 1
      do while (mod(remaining, 2) == 0)
         remaining = remaining / 2
      end do
      do while (mod(remaining, 3) == 0)
         remaining = remaining / 3
      end do
      do while (mod(remaining, 5) == 0)
         remaining = remaining / 5
      end do
      do while (mod(remaining, 7) == 0)
         remaining = remaining / 7
      end do
      !! Rader is attractive when its cyclic convolution length is composed
      !! only of the small radices for which ISFFT has efficient mixed kernels.
      Should_Use_Rader = remaining == 1
   end function Should_Use_Rader

   integer function Modular_Power(base, exponent, modulus) result(value)
      integer, intent(in) :: base, exponent, modulus
      integer :: power
      integer(ISFFT_INT64_KIND) :: factor64, value64, modulus64

      value64 = 1_ISFFT_INT64_KIND
      modulus64 = int(modulus, ISFFT_INT64_KIND)
      factor64 = modulo(int(base, ISFFT_INT64_KIND), modulus64)
      power = exponent
      do while (power > 0)
         if (mod(power, 2) == 1) value64 = modulo(value64 * factor64, modulus64)
         power = power / 2
         if (power > 0) factor64 = modulo(factor64 * factor64, modulus64)
      end do
      value = int(value64)
   end function Modular_Power

   integer function Find_Primitive_Root(prime) result(root)
      integer, intent(in) :: prime
      integer, allocatable :: factors(:)
      integer :: candidate, i, unique_count

      call Factorize_Length(prime - 1, factors)
      unique_count = 0
      do i = 1, size(factors)
         if (i == 1) then
            unique_count = unique_count + 1
            factors(unique_count) = factors(i)
         else if (factors(i) /= factors(i - 1)) then
            unique_count = unique_count + 1
            factors(unique_count) = factors(i)
         end if
      end do

      root = 0
      do candidate = 2, prime - 1
         do i = 1, unique_count
            if (Modular_Power(candidate, (prime - 1) / factors(i), prime) == 1) exit
         end do
         if (i > unique_count) then
            root = candidate
            exit
         end if
      end do
      deallocate (factors)
      if (root == 0) error stop 'Find_Primitive_Root: no primitive root found'
   end function Find_Primitive_Root

   subroutine Build_Rader_Plan(transform_length, slot)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: transform_length
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer, allocatable :: factors(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: temporary_twc(:), temporary_tws(:)
      integer :: q, root, j, k, index_value
      real(kind=ISFFT_REAL_KIND) :: angle, hr, hi, wr, wi
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      if (.not. Is_Prime(transform_length)) return
      q = transform_length - 1
      root = Find_Primitive_Root(transform_length)
      slot%rader_n = q
      slot%rader_is_radix2 = Is_Power_Of_Two(q)

      call Factorize_Length(q, factors)
      slot%rader_num_stages = size(factors)
      allocate (slot%rader_stage_radix(slot%rader_num_stages))
      slot%rader_stage_radix = factors
      deallocate (factors)

      allocate (slot%rader_perm(q), slot%rader_out_perm(q))
      index_value = 1
      do j = 0, q - 1
         slot%rader_perm(j + 1) = index_value
         slot%rader_out_perm(q - j) = index_value
         index_value = &
            int(modulo(int(index_value, ISFFT_INT64_KIND) * int(root, ISFFT_INT64_KIND), &
               int(transform_length, ISFFT_INT64_KIND)))
      end do

      allocate (slot%rader_twc(q), slot%rader_tws(q))
      do k = 0, q - 1
         slot%rader_twc(k + 1) = cos(2.0_ISFFT_REAL_KIND * pi * real(k, ISFFT_REAL_KIND) / real(q, ISFFT_REAL_KIND))
         slot%rader_tws(k + 1) = sin(2.0_ISFFT_REAL_KIND * pi * real(k, ISFFT_REAL_KIND) / real(q, ISFFT_REAL_KIND))
      end do
      if (slot%rader_is_radix2) then
         call Build_Rev_Twiddle(q, slot%rader_rev, temporary_twc, temporary_tws)
         deallocate (temporary_twc, temporary_tws)
      end if

      allocate (slot%rader_kernel_pos_r(q), slot%rader_kernel_pos_i(q))
      allocate (slot%rader_kernel_neg_r(q), slot%rader_kernel_neg_i(q))
      allocate (slot%rader_data_r(q), slot%rader_data_i(q))
      allocate (slot%rader_work_r1(q), slot%rader_work_i1(q))
      allocate (slot%rader_work_r2(q), slot%rader_work_i2(q))
      slot%rader_kernel_pos_r = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_pos_i = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_neg_r = 0.0_ISFFT_REAL_KIND
      slot%rader_kernel_neg_i = 0.0_ISFFT_REAL_KIND
      if (slot%rader_is_radix2) then
         do j = 0, q - 1
            index_value = slot%rader_out_perm(j + 1)
            angle = 2.0_ISFFT_REAL_KIND * pi * real(index_value, ISFFT_REAL_KIND) / &
               real(transform_length, ISFFT_REAL_KIND)
            slot%rader_kernel_pos_r(j + 1) = cos(angle)
            slot%rader_kernel_pos_i(j + 1) = -sin(angle)
            slot%rader_kernel_neg_r(j + 1) = cos(angle)
            slot%rader_kernel_neg_i(j + 1) = sin(angle)
         end do
         call Radix2_Transform_Unscaled(q, slot%rader_kernel_pos_r, slot%rader_kernel_pos_i, +1, slot%rader_rev, &
            slot%rader_twc, slot%rader_tws, slot%rader_work_r1, slot%rader_work_i1)
         call Radix2_Transform_Unscaled(q, slot%rader_kernel_neg_r, slot%rader_kernel_neg_i, +1, slot%rader_rev, &
            slot%rader_twc, slot%rader_tws, slot%rader_work_r1, slot%rader_work_i1)
      else
         do k = 0, q - 1
            do j = 0, q - 1
               index_value = slot%rader_out_perm(j + 1)
               angle = 2.0_ISFFT_REAL_KIND * pi * real(index_value, ISFFT_REAL_KIND) / &
                  real(transform_length, ISFFT_REAL_KIND)
               hr = cos(angle)
               hi = -sin(angle)
               wr = slot%rader_twc(mod(j * k, q) + 1)
               wi = -slot%rader_tws(mod(j * k, q) + 1)
               slot%rader_kernel_pos_r(k + 1) = slot%rader_kernel_pos_r(k + 1) + hr * wr - hi * wi
               slot%rader_kernel_pos_i(k + 1) = slot%rader_kernel_pos_i(k + 1) + hr * wi + hi * wr
               hi = -hi
               slot%rader_kernel_neg_r(k + 1) = slot%rader_kernel_neg_r(k + 1) + hr * wr - hi * wi
               slot%rader_kernel_neg_i(k + 1) = slot%rader_kernel_neg_i(k + 1) + hr * wi + hi * wr
            end do
         end do
      end if
      slot%rader_ready = .true.
   end subroutine Build_Rader_Plan

   subroutine Prepare_Mixed_Factors(transform_length, slot)
      integer, intent(in) :: transform_length
      type(FFT_Precompute_Slot), intent(inout) :: slot
      integer, allocatable :: factors(:)
      integer :: i, current_size
      real(ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      if (allocated(slot%stage_radix)) return
      call Factorize_Length(transform_length, factors)
      slot%num_stages = size(factors)
      allocate (slot%stage_radix(slot%num_stages))
      allocate (slot%stage_size(slot%num_stages))

      current_size = transform_length
      do i = 1, slot%num_stages
         slot%stage_radix(i) = factors(i)
         slot%stage_size(i) = current_size
         current_size = current_size / factors(i)
      end do
      if (.not. allocated(slot%twc)) then
         allocate (slot%twc(transform_length), slot%tws(transform_length))
         do i = 0, transform_length - 1
            slot%twc(i + 1) = &
               cos(2.0_ISFFT_REAL_KIND * pi * real(i, ISFFT_REAL_KIND) / real(transform_length, ISFFT_REAL_KIND))
            slot%tws(i + 1) = &
               sin(2.0_ISFFT_REAL_KIND * pi * real(i, ISFFT_REAL_KIND) / real(transform_length, ISFFT_REAL_KIND))
         end do
      end if
      if (.not. allocated(slot%work_r1)) allocate (slot%work_r1(transform_length), slot%work_i1(transform_length))
      if (.not. allocated(slot%work_r2)) allocate (slot%work_r2(transform_length), slot%work_i2(transform_length))
      if (.not. allocated(slot%work_c1)) allocate (slot%work_c1(transform_length))
      if (.not. allocated(slot%aux_r)) allocate (slot%aux_r(transform_length), slot%aux_i(transform_length))
      deallocate (factors)
   end subroutine ! Prepare_Mixed_Factors

   subroutine Build_Generic_Plan(transform_length, slot)
      integer, intent(in) :: transform_length
      type(FFT_Precompute_Slot), intent(inout) :: slot

      !! ESTIMATE selects the numerical algorithm before building its resources.
      !! Scalar Rader owns its convolution tables; it needs no separate child slot.
      slot%uses_fallback = .true.
      slot%use_rader = Should_Use_Rader(transform_length)
      slot%use_bluestein = .not. slot%use_rader
      if (slot%use_rader) then
         call Build_Rader_Plan(transform_length, slot)
      else
         call Build_Bluestein_Plan(transform_length, slot)
      end if
   end subroutine Build_Generic_Plan

   !! Slot growth can relocate allocatable storage. Reacquire slot views after
   !! recursive preparation before accessing their arrays.
   recursive subroutine Ensure_FFT_Precompute(transform_length, slot_index, real_only)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot_index
      type(FFT_Precompute_Slot), allocatable :: grown_slots(:)
      integer :: capacity, k
      logical, intent(in), optional :: real_only
      logical :: coefficients_only
      real(ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      if (transform_length < 1) error stop 'Ensure_FFT_Precompute: len must be positive'

      coefficients_only = .false.
      if (present(real_only)) coefficients_only = real_only
      slot_index = Find_FFT_Precompute_Slot(transform_length)
      if (slot_index > 0) then
         fft_precompute_use_clock = fft_precompute_use_clock + 1
         fft_precompute_slots(slot_index)%last_used = fft_precompute_use_clock
         if (coefficients_only) then
            if (allocated(fft_precompute_slots(slot_index)%twc) .and. &
               allocated(fft_precompute_slots(slot_index)%work_r1)) &
               return
         else
            if (fft_precompute_slots(slot_index)%complex_ready) return
         end if
      else

         slot_index = fft_precompute_count + 1
         if (.not. allocated(fft_precompute_slots)) then
            allocate (fft_precompute_slots(slot_index))
         else if (slot_index > size(fft_precompute_slots)) then
            capacity = slot_index
            capacity = max(slot_index, 2 * size(fft_precompute_slots))
            allocate (grown_slots(capacity))
            if (fft_precompute_count > 0) grown_slots(1:fft_precompute_count) = &
               fft_precompute_slots(1:fft_precompute_count)
            call move_alloc(grown_slots, fft_precompute_slots)
         end if
         fft_precompute_count = slot_index

         call Reset_FFT_Precompute_Slot(fft_precompute_slots(slot_index))
      end if
      if (.not. allocated(fft_precompute_slots(slot_index)%r2c_re)) &
         allocate (fft_precompute_slots(slot_index)%r2c_re(transform_length), &
         fft_precompute_slots(slot_index)%r2c_im(transform_length))
      fft_precompute_slots(slot_index)%is_radix2 = Is_Power_Of_Two(transform_length)
      if (coefficients_only) then
         !! Even real FFT: N-point reconstruction coefficients, only N/2-point FFT.
         if (.not. allocated(fft_precompute_slots(slot_index)%twc)) then
            allocate (fft_precompute_slots(slot_index)%twc(transform_length), &
               fft_precompute_slots(slot_index)%tws(transform_length))
            do k = 0, transform_length - 1
               fft_precompute_slots(slot_index)%twc(k + 1) = cos(2 * pi * k / transform_length)
               fft_precompute_slots(slot_index)%tws(k + 1) = sin(2 * pi * k / transform_length)
            end do
         end if
         if (.not. allocated(fft_precompute_slots(slot_index)%work_r1)) &
            allocate (fft_precompute_slots(slot_index)%work_r1(transform_length), &
            fft_precompute_slots(slot_index)%work_i1(transform_length))
         if (.not. allocated(fft_precompute_slots(slot_index)%work_c1)) &
            allocate (fft_precompute_slots(slot_index)%work_c1(transform_length))
      else
         if (fft_precompute_slots(slot_index)%is_radix2) then
            if (allocated(fft_precompute_slots(slot_index)%twc)) then
               !! Promote real reconstruction metadata without recomputing its factors.
               call Build_Bit_Reverse(transform_length, fft_precompute_slots(slot_index)%rev)
            else
               call Build_Rev_Twiddle(transform_length, fft_precompute_slots(slot_index)%rev, &
                  fft_precompute_slots(slot_index)%twc, &
                  fft_precompute_slots(slot_index)%tws)
            end if
            allocate (fft_precompute_slots(slot_index)%twiddle_forward(transform_length))
            allocate (fft_precompute_slots(slot_index)%twiddle_backward(transform_length))
            if (transform_length == 1) then
               !! The one-point FFT is identity; there is no negative half-table.
               fft_precompute_slots(slot_index)%twiddle_forward = &
                  cmplx(1.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
               fft_precompute_slots(slot_index)%twiddle_backward = fft_precompute_slots(slot_index)%twiddle_forward
            else
               fft_precompute_slots(slot_index)%twiddle_forward(1:transform_length / 2) = &
                  cmplx(fft_precompute_slots(slot_index)%twc(1:transform_length / 2), - &
                  fft_precompute_slots(slot_index)%tws(1:transform_length / 2), &
                  kind=ISFFT_REAL_KIND)
               fft_precompute_slots(slot_index)%twiddle_backward(1:transform_length / 2) = &
                  cmplx(fft_precompute_slots(slot_index)%twc(1:transform_length / 2), &
                  fft_precompute_slots(slot_index)%tws(1:transform_length / 2), &
                  kind=ISFFT_REAL_KIND)
               fft_precompute_slots(slot_index)%twiddle_forward(transform_length / 2 + 1:transform_length) = - &
                  fft_precompute_slots(slot_index)%twiddle_forward(1:transform_length / 2)
               fft_precompute_slots(slot_index)%twiddle_backward(transform_length / 2 + 1:transform_length) = - &
                  fft_precompute_slots(slot_index)%twiddle_backward(1:transform_length / 2)
            end if
            if (.not. allocated(fft_precompute_slots(slot_index)%work_r1)) &
               allocate (fft_precompute_slots(slot_index)%work_r1(transform_length))
            if (.not. allocated(fft_precompute_slots(slot_index)%work_i1)) &
               allocate (fft_precompute_slots(slot_index)%work_i1(transform_length))
            allocate (fft_precompute_slots(slot_index)%work_r2(transform_length))
            allocate (fft_precompute_slots(slot_index)%work_i2(transform_length))
            if (.not. allocated(fft_precompute_slots(slot_index)%work_c1)) &
               allocate (fft_precompute_slots(slot_index)%work_c1(transform_length))
            allocate (fft_precompute_slots(slot_index)%aux_r(transform_length))
            allocate (fft_precompute_slots(slot_index)%aux_i(transform_length))
         else
            if (.not. allocated(fft_precompute_slots(slot_index)%work_c1)) &
               allocate (fft_precompute_slots(slot_index)%work_c1(transform_length))
            allocate (fft_precompute_slots(slot_index)%aux_r(transform_length), &
               fft_precompute_slots(slot_index)%aux_i(transform_length))
            call Build_Generic_Plan(transform_length, fft_precompute_slots(slot_index))
         end if
         fft_precompute_slots(slot_index)%complex_ready = .true.
      end if
      fft_precompute_slots(slot_index)%transform_length = transform_length
      fft_precompute_slots(slot_index)%active = .true.
      call Reserve_FFT_Base(slot_index)
      fft_precompute_use_clock = fft_precompute_use_clock + 1
      fft_precompute_slots(slot_index)%last_used = fft_precompute_use_clock
      !! Prepare subordinate transforms before any numerical kernel is entered.
      block
      end block
   end subroutine Ensure_FFT_Precompute

   subroutine Ensure_Rader_Precompute(transform_length, slot_index)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot_index

      call Ensure_FFT_Precompute(transform_length, slot_index)
      if (.not. Is_Prime(transform_length)) error stop 'Ensure_Rader_Precompute: length must be prime'
      if (.not. fft_precompute_slots(slot_index)%rader_ready) then
         call Build_Rader_Plan(transform_length, fft_precompute_slots(slot_index))
      end if
      block
         integer :: child
         if (transform_length > 2) call Ensure_FFT_Precompute(transform_length - 1, child)
      end block
   end subroutine Ensure_Rader_Precompute

   subroutine Ensure_R2C_Precompute(transform_length, slot)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot
      integer :: child

      call Ensure_FFT_Precompute(transform_length, slot, real_only=mod(transform_length, 2) == 0)
      if (mod(transform_length, 2) == 0) call Ensure_FFT_Precompute(transform_length / 2, child)
   end subroutine ! Ensure_R2C_Precompute

   subroutine Ensure_Real_Transform_Precompute(packed_len, slot_index)
      use ISFFT_Precision, only: ISFFT_REAL_KIND
      integer, intent(in) :: packed_len
      integer, intent(out), optional :: slot_index
      integer :: resolved_slot, k
      real(kind=ISFFT_REAL_KIND), parameter :: pi = acos(-1.0_ISFFT_REAL_KIND)

      if (mod(packed_len, 2) == 0) then
         call Ensure_R2C_Precompute(packed_len, resolved_slot)
      else
         call Ensure_FFT_Precompute(packed_len, resolved_slot)
      end if
      if (.not. allocated(fft_precompute_slots(resolved_slot)%real_twc)) then
         allocate (fft_precompute_slots(resolved_slot)%real_twc(packed_len + 1))
         allocate (fft_precompute_slots(resolved_slot)%real_tws(packed_len + 1))
         do k = 0, packed_len
            fft_precompute_slots(resolved_slot)%real_twc(k + 1) = &
               cos(pi * real(k, ISFFT_REAL_KIND) / real(packed_len, ISFFT_REAL_KIND))
            fft_precompute_slots(resolved_slot)%real_tws(k + 1) = &
               sin(pi * real(k, ISFFT_REAL_KIND) / real(packed_len, ISFFT_REAL_KIND))
         end do
      end if
      if (present(slot_index)) slot_index = resolved_slot
   end subroutine Ensure_Real_Transform_Precompute

   subroutine Add_Unique_Length(domain_lengths, count, transform_length)
      integer, intent(inout) :: domain_lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: transform_length

      if (transform_length < 1) return
      if (count > 0) then
         if (any(domain_lengths(1:count) == transform_length)) return
      end if
      if (count >= size(domain_lengths)) error stop 'FFT_Initialize: precompute capacity is too small'
      count = count + 1
      domain_lengths(count) = transform_length
   end subroutine Add_Unique_Length

   subroutine Add_Periodic_Lengths(lengths, count, grid_len)
      integer, intent(inout) :: lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: grid_len

      call Add_Unique_Length(lengths, count, grid_len)
      if (mod(grid_len, 2) == 0) call Add_Unique_Length(lengths, count, grid_len / 2)
   end subroutine Add_Periodic_Lengths

   subroutine Add_Nonperiodic_Lengths(lengths, count, grid_len)
      integer, intent(inout) :: lengths(:)
      integer, intent(inout) :: count
      integer, intent(in) :: grid_len

      if (grid_len < 2) error stop 'FFT_Initialize: nonperiodic grid length must be at least two'
      call Add_Unique_Length(lengths, count, grid_len - 1)
   end subroutine Add_Nonperiodic_Lengths

   subroutine FFT_Initialize(nx, ny, nz, boundary_types, active_dimension)
      integer, intent(in) :: nx, ny, nz
      integer, intent(in), optional :: boundary_types(3)
      integer, intent(in), optional :: active_dimension
      integer :: grid_lengths(3)
      integer, allocatable :: init_lengths(:)
      integer :: i, slot_idx, init_count, dimension_count

      grid_lengths = [nx, ny, nz]
      init_count = 0
      dimension_count = 3
      if (present(active_dimension)) dimension_count = active_dimension
      if (dimension_count < 1 .or. dimension_count > 3) then
         error stop 'FFT_Initialize: active_dimension must be 1, 2, or 3'
      end if

      call Clear_FFT_Precompute_Cache()
      if (allocated(real_workspaces)) deallocate (real_workspaces)
      if (allocated(prepared_kernel_plans)) deallocate (prepared_kernel_plans)
      !! A dimensions can contribute n, n/2 (even periodic R2C), and n-1
      !! (nonperiodic packed transform); this is an exact semantic bound.
      allocate (init_lengths(3 * dimension_count))
      init_lengths = 0

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
      !! Allocate the exact known semantic set before constructing any slot.
      !! This avoids moving plans that already own large allocatable workspaces.
      allocate (fft_precompute_slots(init_count))
      do i = 1, init_count
         call Ensure_FFT_Precompute(init_lengths(i), slot_idx)
         call Reserve_R2C_Workspace(init_lengths(i), 128)
         call Reserve_FFT_Axis(init_lengths(i))
      end do
      deallocate (init_lengths)

      do i = 1, dimension_count
         if (present(boundary_types)) then
            if (boundary_types(i) == 1 .or. boundary_types(i) == 2) then
               call Ensure_Real_Transform_Precompute(grid_lengths(i) - 1)
            end if
         else
            call Ensure_Real_Transform_Precompute(grid_lengths(i) - 1)
         end if
      end do

   end subroutine FFT_Initialize

   subroutine FFT_Finalize()

      call Clear_FFT_Precompute_Cache()
      if (allocated(real_workspaces)) deallocate (real_workspaces)
      if (allocated(prepared_kernel_plans)) deallocate (prepared_kernel_plans)
   end subroutine FFT_Finalize

   subroutine Prepare_Real_Workspace(transform_length, line_count)
      integer, intent(in) :: transform_length, line_count
      type(Real_Transform_Workspace), allocatable :: grown(:)
      integer :: i, count, q, p, fft_len
      logical :: build_indices

      if (mod(transform_length, 2) == 0) call Reserve_R2C_Workspace(transform_length, line_count)
      if (.not. allocated(real_workspaces)) allocate (real_workspaces(0))
      do i = 1, size(real_workspaces)
         if (real_workspaces(i)%transform_length == transform_length) exit
      end do
      if (i > size(real_workspaces)) then
         count = size(real_workspaces)
         allocate (grown(count + 1))
         grown(1:count) = real_workspaces
         call move_alloc(grown, real_workspaces)
      end if
      associate (work => real_workspaces(i))
         if (work%line_count >= line_count .and. work%transform_length == transform_length) return
         build_indices = .not. allocated(work%real_work)
         if (allocated(work%real_work)) then
            deallocate (work%real_work, work%spectrum, work%odd_packed, work%batch_storage)
         else
            allocate (work%dst_pack_real_index(transform_length), work%dst_pack_imag_index(transform_length), &
               work%dst_pack_real_sign(transform_length), &
               work%dst_pack_imag_sign(transform_length), work%dct_pack_real_index(transform_length), &
               work%dct_pack_imag_index(transform_length))
         end if
         work%transform_length = transform_length
         work%line_count = max(1, line_count)
         allocate (work%real_work(transform_length * work%line_count), &
            work%spectrum((transform_length / 2 + 1) * work%line_count), &
            work%odd_packed(transform_length * ((work%line_count + 1) / 2)), &
            work%batch_storage(transform_length * work%line_count))
         if (build_indices) then
            fft_len = 2 * transform_length
            do q = 0, transform_length - 1
               p = 2 * q
               work%dct_pack_real_index(q + 1) = min(p, fft_len - p) + 1
               work%dst_pack_real_index(q + 1) = max(1, min(p, fft_len - p))
               work%dst_pack_real_sign(q + 1) = 1
               if (p > transform_length) work%dst_pack_real_sign(q + 1) = -1
               if (p == 0 .or. p == transform_length) then
                  work%dst_pack_real_sign(q + 1) = 0
                  work%dst_pack_real_index(q + 1) = 1
               end if
               p = 2 * q + 1
               work%dct_pack_imag_index(q + 1) = min(p, fft_len - p) + 1
               work%dst_pack_imag_index(q + 1) = max(1, min(p, fft_len - p))
               work%dst_pack_imag_sign(q + 1) = 1
               if (p > transform_length) work%dst_pack_imag_sign(q + 1) = -1
               if (p == transform_length) then
                  work%dst_pack_imag_sign(q + 1) = 0
                  work%dst_pack_imag_index(q + 1) = 1
               end if
            end do
         end if
      end associate
   end subroutine ! Prepare_Real_Workspace

   subroutine Prepare_Complex_Candidates(transform_length, slot)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot
      integer :: child, q

      call Ensure_FFT_Precompute(transform_length, slot)
      if (fft_precompute_slots(slot)%candidates_ready) return
      if (fft_precompute_slots(slot)%is_radix2) then
         fft_precompute_slots(slot)%is_radix4 = Is_Power_Of_Four(transform_length)
         if (fft_precompute_slots(slot)%is_radix4) &
            call Build_Radix4_Reverse(transform_length, fft_precompute_slots(slot)%rev_radix4)
      else
         call Prepare_Mixed_Factors(transform_length, fft_precompute_slots(slot))
         if (fft_precompute_slots(slot)%bluestein_m == 0) &
            call Build_Bluestein_Plan(transform_length, fft_precompute_slots(slot))
      end if
      fft_precompute_slots(slot)%batch_ready = .true.
      call Reserve_FFT_Base(slot)
      if (ISFFT_Is_Three_Times_Power_Of_Two(transform_length)) call Ensure_FFT_Precompute(transform_length / 3, child)
      if (fft_precompute_slots(slot)%rader_ready) then
         q = transform_length - 1
         call Ensure_FFT_Precompute(q, child)
         if (.not. fft_precompute_slots(child)%is_radix2) call Prepare_Mixed_Factors(q, fft_precompute_slots(child))
         fft_precompute_slots(child)%batch_ready = .true.
         call Reserve_FFT_Base(child)
      end if
      fft_precompute_slots(slot)%candidates_ready = .true.
   end subroutine ! Prepare_Complex_Candidates

   subroutine Prepare_R2C_Candidates(transform_length, line_count, slot)
      integer, intent(in) :: transform_length, line_count
      integer, intent(out) :: slot
      integer :: child, third

      call Ensure_R2C_Precompute(transform_length, slot)
      if (mod(transform_length, 2) == 0) then
         call Ensure_FFT_Precompute(transform_length / 2, child)
         !! Generic real candidate transforms N complex values; power-of-two real
         !! transforms need only reconstruction metadata at N, not an N-point FFT.
         if (.not. Is_Power_Of_Two(transform_length)) then
            call Ensure_FFT_Precompute(transform_length, slot)
            fft_precompute_slots(slot)%batch_ready = .true.
            call Reserve_FFT_Base(slot)
            if (ISFFT_Is_Three_Times_Power_Of_Two(transform_length)) then
               call Prepare_Mixed_Factors(transform_length, fft_precompute_slots(slot))
               call Ensure_FFT_Precompute(transform_length / 3, third)
            end if
         end if
         if (ISFFT_Is_Three_Times_Power_Of_Two(transform_length / 2)) then
            call Prepare_Mixed_Factors(transform_length / 2, fft_precompute_slots(child))
            fft_precompute_slots(child)%batch_ready = .true.
            call Reserve_FFT_Base(child)
            call Ensure_FFT_Precompute(transform_length / 6, third)
         end if
      end if
      call Reserve_R2C_Workspace(transform_length, line_count)
   end subroutine ! Prepare_R2C_Candidates

   subroutine Build_Poisson_Operator(op, spectral_layout, global_shape, boundary_types, dimensions, domain_lengths, &
      discretization_scheme)
      type(Poisson_Operator), intent(inout) :: op
      type(Layout), intent(in) :: spectral_layout
      integer, intent(in) :: global_shape(3), boundary_types(3), dimensions, discretization_scheme
      real(kind=ISFFT_REAL_KIND), intent(in) :: domain_lengths(3)
      type(ISFFT_Scheme_Coefficients) :: coeff
      integer :: rank, ierr, n(3), off(3), sh(3), i, j, k, a, modes(3), mode
      real(kind=ISFFT_REAL_KIND) :: lambda
      real(kind=ISFFT_REAL_KIND), allocatable :: axis_lambda(:, :)
      logical :: boundary

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call Layout_Box(spectral_layout, rank, n, off)
      sh = n(spectral_layout%axis_order)
      call Destroy_Poisson_Operator(op)
      allocate (op%multiplier(sh(1), sh(2), sh(3)))
      call ISFFT_Get_Scheme_Coefficients(discretization_scheme, coeff)
      !! The Laplacian symbol is separable. Evaluate each local 1D mode once,
      !! then build the same full multiplier used by the unchanged apply loop.
      allocate (axis_lambda(maxval(n), dimensions))
      do a = 1, dimensions
         do i = 1, n(a)
            mode = off(a) + i - 1
            if (boundary_types(a) == 0 .and. 2 * mode >= global_shape(a)) mode = mode - global_shape(a)
            axis_lambda(i, a) = &
               ISFFT_Lambda_For_BC(boundary_types(a), mode, global_shape(a), domain_lengths(a), &
               discretization_scheme, coeff)
         end do
      end do
      do k = 1, sh(3)
         do j = 1, sh(2)
            do i = 1, sh(1)
               modes(spectral_layout%axis_order) = [i - 1, j - 1, k - 1]
               modes = modes + off
               lambda = 0
               boundary = .false.
               do a = 1, dimensions
                  mode = modes(a)
                  if (boundary_types(a) == 1 .and. (mode == 0 .or. mode == global_shape(a) - 1)) boundary = .true.
                  lambda = lambda + axis_lambda(mode - off(a) + 1, a)
               end do
               op%multiplier(i, j, k) = 0
               if (.not. boundary .and. lambda /= 0) op%multiplier(i, j, k) = -1.0_ISFFT_REAL_KIND / lambda
            end do
         end do
      end do
   end subroutine ! Build_Poisson_Operator

end module ISFFT_Precompute
