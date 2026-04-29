module ISFFT_FFT_Pre
   use, intrinsic :: iso_fortran_env, only: output_unit
   use ISFFT_Buffer, only: &
      FFT_PRECOMPUTE_MAX, FFT_Precompute_Slot, &
      fft_precompute_slots, fft_precompute_count, fft_precompute_next_slot
   implicit none

   private
   public :: FFT_Initialize, FFT_Finalize, Ensure_FFT_Precompute, &
      Ensure_Generic_N, Ensure_Generic_M, Build_Rev_Twiddle, &
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
      if (allocated(slot%twc)) deallocate(slot%twc)
      if (allocated(slot%tws)) deallocate(slot%tws)

      slot%len = 0
      slot%active = .false.
   end subroutine Reset_FFT_Precompute_Slot
   !----------------------------------------------------------------------
   subroutine Clear_FFT_Precompute_Cache()
      integer :: i

      do i = 1, FFT_PRECOMPUTE_MAX
         call Reset_FFT_Precompute_Slot(fft_precompute_slots(i))
      end do

      fft_precompute_count = 0
      fft_precompute_next_slot = 1
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
   subroutine Ensure_FFT_Precompute(len, slot_idx)
      integer, intent(in) :: len
      integer, intent(out) :: slot_idx

      call Validate_FFT_Length(len, 'Ensure_FFT_Precompute')

      slot_idx = Find_FFT_Precompute_Slot(len)
      if (slot_idx > 0) return

      if (fft_precompute_count < FFT_PRECOMPUTE_MAX) then
         slot_idx = fft_precompute_count + 1
         fft_precompute_count = slot_idx
         if (fft_precompute_count < FFT_PRECOMPUTE_MAX) then
            fft_precompute_next_slot = fft_precompute_count + 1
         else
            fft_precompute_next_slot = 1
         end if
      else
         slot_idx = fft_precompute_next_slot
         fft_precompute_next_slot = mod(slot_idx, FFT_PRECOMPUTE_MAX) + 1
      end if

      call Build_Rev_Twiddle(len, fft_precompute_slots(slot_idx)%rev, &
         fft_precompute_slots(slot_idx)%twc, fft_precompute_slots(slot_idx)%tws)
      fft_precompute_slots(slot_idx)%len = len
      fft_precompute_slots(slot_idx)%active = .true.
      ! call FFT_Print_Precompute_Cache()
   end subroutine Ensure_FFT_Precompute
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
   subroutine FFT_Initialize(nx, ny, nz)
      integer, intent(in) :: nx, ny, nz
      integer :: init_lengths(6)
      integer :: i, slot_idx

      init_lengths = [nx, ny, nz, nx/2, ny/2, nz/2]

      call Clear_FFT_Precompute_Cache()

      do i = 1, size(init_lengths)
         call Ensure_FFT_Precompute(init_lengths(i), slot_idx)
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

         write(out_unit,'(A,I0,A,A,A,I0,1X,A)') 'slot ', i, ': ', trim(slot_state), &
            ', len=', fft_precompute_slots(i)%len, trim(next_marker)
      end do
   end subroutine FFT_Print_Precompute_Cache
   !----------------------------------------------------------------------

end module ISFFT_FFT_Pre
