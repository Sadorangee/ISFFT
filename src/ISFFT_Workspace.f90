!------------------------------------------------------------------------------
! ISFFT_Workspace
! Numerical tables, cached records and reusable field and communication storage.
!------------------------------------------------------------------------------
module ISFFT_Workspace
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND
   implicit none
   private
   public :: FFT_Precompute_Slot, Real_Transform_Workspace, Field_Buffer, Spectral_Field, Poisson_Operator
   public :: Access_Candidate, Sample_Scores

   public :: fft_precompute_slots, fft_precompute_count, fft_precompute_use_clock, real_workspaces, sample_cache

   public :: poisson_field, redistribution_field, poisson_symbol

   public :: axis_scratch, axis_real, axis_imag, block_send_r, block_recv_r
   public :: block_tile_r, block_send_z, block_recv_z, block_tile_z, block_counts
   public :: block_displs, halo_work, workspace_growths

   public :: send_real_buffer, recv_real_buffer, send_complex_buffer, recv_complex_buffer, workspace_sealed

   public :: Clear_Field, Park_Field, Get_Real, Get_Complex, Load_Field
   public :: Destroy_Poisson_Operator

   public :: Reserve_Real, Reserve_Complex, Reserve_R, Reserve_Z, Release_Execution_Scratch

   public :: Require_FFT_Precompute, Require_Rader_Precompute, Require_Real_Transform_Precompute, Require_Real_Workspace

   public :: Seal_Layout_Workspace, Release_Layout_Workspace, Bound_Layout_Workspace

   !! FFT precompute records
   type :: FFT_Precompute_Slot
      integer :: transform_length = 0
      integer :: last_used = 0
      logical :: active = .false.
      logical :: complex_ready = .false.
      logical :: batch_ready = .false., candidates_ready = .false.
      logical :: is_radix2 = .false.
      logical :: is_radix4 = .false.
      logical :: uses_fallback = .false.
      logical :: use_bluestein = .false.
      logical :: use_rader = .false.
      logical :: rader_ready = .false.
      logical :: rader_is_radix2 = .false.
      integer :: num_stages = 0
      integer :: bluestein_m = 0
      integer :: rader_n = 0
      integer :: rader_num_stages = 0
      integer :: complex_strategy_count = 0
      integer, allocatable :: complex_strategy_lines(:)
      integer, allocatable :: complex_strategy_method(:)
      integer, allocatable :: complex_strategy_width(:)
      integer :: r2c_strategy_count = 0
      integer, allocatable :: r2c_strategy_lines(:)
      integer, allocatable :: r2c_strategy_method(:)
      integer, allocatable :: r2c_strategy_width(:)
      integer, allocatable :: rev(:)
      integer, allocatable :: rev_radix4(:)
      integer, allocatable :: stage_radix(:), stage_size(:)
      integer, allocatable :: rader_perm(:), rader_out_perm(:), rader_stage_radix(:), rader_rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: twc(:), tws(:)
      complex(kind=ISFFT_REAL_KIND), allocatable :: twiddle_forward(:), twiddle_backward(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: real_twc(:), real_tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: work_r1(:), work_i1(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: work_r2(:), work_i2(:)
      complex(kind=ISFFT_REAL_KIND), allocatable :: work_c1(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: aux_r(:), aux_i(:), r2c_re(:), r2c_im(:)
      integer, allocatable :: bluestein_rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_twc(:), bluestein_tws(:)
      complex(kind=ISFFT_REAL_KIND), allocatable :: bluestein_twiddle_forward(:), bluestein_twiddle_backward(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: chirp_c(:), chirp_s(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: kernel_pos_r(:), kernel_pos_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: kernel_neg_r(:), kernel_neg_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_r(:), bluestein_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_work_r(:), bluestein_work_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_twc(:), rader_tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_kernel_pos_r(:), rader_kernel_pos_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_kernel_neg_r(:), rader_kernel_neg_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_data_r(:), rader_data_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_work_r1(:), rader_work_i1(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_work_r2(:), rader_work_i2(:)
      complex(kind=ISFFT_REAL_KIND), allocatable :: line_work(:), split_work(:), batch_work(:), mixed_1(:), &
         mixed_2(:), &
         mixed_3(:), rader_batch(:), axis_batch(:), axis_lines(:), r2c_pack(:), r2c_natural(:), r2c_generic(:), &
         r2c_lines(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: scalar_re(:), scalar_im(:)
   end type FFT_Precompute_Slot

   !! R2R scratch and extension maps are shared sequentially by DST and DCT.
   !! Keyed by interval length; each entry reserves the maximum requested line count.
   type :: Real_Transform_Workspace
      integer :: transform_length = 0, line_count = 0
      real(kind=ISFFT_REAL_KIND), allocatable :: real_work(:)
      complex(kind=ISFFT_REAL_KIND), allocatable :: spectrum(:), odd_packed(:), batch_storage(:)
      integer, allocatable :: dst_pack_real_index(:), dst_pack_imag_index(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: dst_pack_real_sign(:), dst_pack_imag_sign(:)
      integer, allocatable :: dct_pack_real_index(:), dct_pack_imag_index(:)
   end type
   type(Real_Transform_Workspace), allocatable, target, save :: real_workspaces(:)
   type(FFT_Precompute_Slot), allocatable, target, save :: fft_precompute_slots(:)
   integer, save :: fft_precompute_count = 0
   integer, save :: fft_precompute_use_clock = 0

   !! Shared execution scratch
   !! Shared scratch assumes sequential execution of library calls within each rank.
   ! Execution scratch is shared sequentially, reserved during preparation only.
   complex(kind=ISFFT_REAL_KIND), allocatable, target, save :: axis_scratch(:)
   real(kind=ISFFT_REAL_KIND), allocatable, target, save :: axis_real(:), axis_imag(:)
   real(kind=ISFFT_REAL_KIND), allocatable, target, save :: block_send_r(:), block_recv_r(:), block_tile_r(:)
   complex(kind=ISFFT_REAL_KIND), allocatable, target, save :: block_send_z(:), block_recv_z(:), block_tile_z(:)
   integer, allocatable, save :: block_counts(:), block_displs(:)
   real(kind=ISFFT_REAL_KIND), allocatable, target, save :: halo_work(:, :)
   integer, save :: workspace_growths = 0

   !! Field storage and allocation ownership
   type :: Field_Buffer
      real(kind=ISFFT_REAL_KIND), allocatable :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: z(:, :, :)
      type(Field_Buffer), pointer :: next => null()
   end type

   type :: Spectral_Field
      logical :: prepared = .false.
      real(kind=ISFFT_REAL_KIND), allocatable :: r(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: z(:, :, :)
      type(Field_Buffer), pointer :: pool => null()
   end type

   type :: Access_Candidate
      logical :: direct_access = .false.
      integer :: batch_width = 1
   end type

   type Sample_Scores
      integer :: transform_length, sample_line_count, candidate_count, planning_level, warmup_count, repeat_count
      integer, allocatable :: methods(:), widths(:)
      real(ISFFT_REAL_KIND), allocatable :: scores(:)
   end type

   type :: Poisson_Operator
      real(kind=ISFFT_REAL_KIND), allocatable :: multiplier(:, :, :)
   end type
   type(Spectral_Field), save :: redistribution_field
   type(Spectral_Field), target, save :: poisson_field
   type(Poisson_Operator), target, save :: poisson_symbol
   type(Sample_Scores), allocatable :: sample_cache(:)

   !! Shared transport buffers
   real(kind=ISFFT_REAL_KIND), allocatable, save, asynchronous :: send_real_buffer(:), recv_real_buffer(:)
   complex(kind=ISFFT_REAL_KIND), allocatable, save, asynchronous :: send_complex_buffer(:), recv_complex_buffer(:)
   logical, save :: workspace_sealed = .false.

contains

   !------------------------------------------------------------------------------
   !! Storage capacity management and prepared-data checks
   !------------------------------------------------------------------------------
   subroutine Reserve_Real(storage, required_capacity)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: storage(:)
      integer(ISFFT_INT64_KIND), intent(in) :: required_capacity

      if (required_capacity <= 0) return
      if (allocated(storage)) then
         if (size(storage, kind=ISFFT_INT64_KIND) >= required_capacity) return
         deallocate (storage)
      end if
      allocate (storage(required_capacity))
      workspace_growths = workspace_growths + 1
   end subroutine ! Reserve_Real

   subroutine Reserve_Complex(storage, required_capacity)
      complex(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: storage(:)
      integer(ISFFT_INT64_KIND), intent(in) :: required_capacity

      if (required_capacity <= 0) return
      if (allocated(storage)) then
         if (size(storage, kind=ISFFT_INT64_KIND) >= required_capacity) return
         deallocate (storage)
      end if
      allocate (storage(required_capacity))
      workspace_growths = workspace_growths + 1
   end subroutine ! Reserve_Complex

   subroutine Release_Execution_Scratch()

      if (allocated(halo_work)) deallocate (halo_work)
      if (allocated(axis_scratch)) deallocate (axis_scratch)
      if (allocated(axis_real)) deallocate (axis_real, axis_imag)
      if (allocated(block_send_r)) deallocate (block_send_r, block_recv_r, block_tile_r)
      if (allocated(block_send_z)) deallocate (block_send_z, block_recv_z, block_tile_z)
      if (allocated(block_counts)) deallocate (block_counts, block_displs)
   end subroutine ! Release_Execution_Scratch

   !! Execution requires prepared factors and scratch. Keep allocation in preparation.
   subroutine Require_FFT_Precompute(transform_length, slot)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot
      integer :: i

      if (allocated(fft_precompute_slots)) then
         do i = 1, size(fft_precompute_slots)
            if (fft_precompute_slots(i)%active .and. fft_precompute_slots(i)%transform_length == transform_length) then
               slot = i
               return
            end if
         end do
      end if
      write (*, *) 'ISFFT: missing FFT preparation for length ', transform_length
      error stop 'Call ISFFT_Initialize or ISFFT_Prepare_Kernel before executing a kernel'
   end subroutine ! Require_FFT_Precompute

   subroutine Require_Rader_Precompute(transform_length, slot)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot

      call Require_FFT_Precompute(transform_length, slot)
      if (.not. fft_precompute_slots(slot)%rader_ready) &
         error stop 'ISFFT: Rader data missing; initialize with ISFFT_Prepare_Kernel first'
   end subroutine ! Require_Rader_Precompute

   subroutine Require_Real_Transform_Precompute(transform_length, slot)
      integer, intent(in) :: transform_length
      integer, intent(out) :: slot

      call Require_FFT_Precompute(transform_length, slot)
      if (.not. allocated(fft_precompute_slots(slot)%real_twc)) &
         error stop 'ISFFT: real-transform factors missing; call ISFFT_Prepare_Kernel first'
   end subroutine ! Require_Real_Transform_Precompute

   subroutine Require_Real_Workspace(transform_length, line_count, work)
      integer, intent(in) :: transform_length, line_count
      type(Real_Transform_Workspace), pointer, intent(out) :: work
      integer :: i

      if (allocated(real_workspaces)) then
         do i = 1, size(real_workspaces)
            if (real_workspaces(i)%transform_length /= transform_length) cycle
            if (real_workspaces(i)%line_count < line_count) exit
            work => real_workspaces(i)
            return
         end do
      end if
      write (*, *) 'ISFFT: missing R2R workspace for interval/lines ', transform_length, line_count
      error stop 'Call ISFFT_Prepare_Kernel with sufficient line capacity before execution'
   end subroutine ! Require_Real_Workspace

   !------------------------------------------------------------------------------
   !! Field pool lifetime and execution-time reuse
   !------------------------------------------------------------------------------
   subroutine Clear_Field(field)
      type(Spectral_Field), intent(inout) :: field
      type(Field_Buffer), pointer :: node, next

      if (allocated(field%r)) deallocate (field%r)
      if (allocated(field%z)) deallocate (field%z)
      node => field%pool
      do while (associated(node))
         next => node%next
         deallocate (node)
         node => next
      end do
      nullify (field%pool)
      field%prepared = .false.
   end subroutine ! Clear_Field

   !! move_alloc transfers ownership into the pool. Reacquire views before reuse.
   subroutine Park_Field(field)
      type(Spectral_Field), intent(inout) :: field
      type(Field_Buffer), pointer :: node

      if (.not. allocated(field%r) .and. .not. allocated(field%z)) return
      node => field%pool
      do while (associated(node))
         if (.not. allocated(node%r) .and. .not. allocated(node%z)) exit
         node => node%next
      end do
      if (.not. associated(node)) then
         if (field%prepared) error stop 'Park_Field: missing reserved pool node'
         allocate (node)
         node%next => field%pool
         field%pool => node
      end if
      if (allocated(field%r)) call move_alloc(field%r, node%r)
      if (allocated(field%z)) call move_alloc(field%z, node%z)
   end subroutine ! Park_Field

   subroutine Get_Real(field, local_shape, r, lower)
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: local_shape(3)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(out) :: r(:, :, :)
      integer, intent(in), optional :: lower
      integer :: lo
      type(Field_Buffer), pointer :: node

      lo = 1
      if (present(lower)) lo = lower
      node => field%pool
      do while (associated(node))
         if (allocated(node%r)) then
            if (all(shape(node%r) == local_shape) .and. all(lbound(node%r) == lo)) then
               call move_alloc(node%r, r)
               return
            end if
         end if
         node => node%next
      end do
      if (field%prepared) error stop 'unreserved real execution field'
      allocate (r(lo:lo + local_shape(1) - 1, lo:lo + local_shape(2) - 1, lo:lo + local_shape(3) - 1))
      workspace_growths = workspace_growths + 1
   end subroutine ! Get_Real

   subroutine Get_Complex(field, local_shape, z)
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: local_shape(3)
      complex(kind=ISFFT_REAL_KIND), allocatable, intent(out) :: z(:, :, :)
      type(Field_Buffer), pointer :: node

      node => field%pool
      do while (associated(node))
         if (allocated(node%z)) then
            if (all(shape(node%z) == local_shape)) then
               call move_alloc(node%z, z)
               return
            end if
         end if
         node => node%next
      end do
      if (field%prepared) error stop 'unreserved complex execution field'
      allocate (z(local_shape(1), local_shape(2), local_shape(3)))
      workspace_growths = workspace_growths + 1
   end subroutine ! Get_Complex

   !------------------------------------------------------------------------------
   !! Physical input loading and layout movement
   !------------------------------------------------------------------------------
   subroutine Load_Field(field, input)
      type(Spectral_Field), intent(inout) :: field
      real(kind=ISFFT_REAL_KIND), intent(in) :: input(:, :, :)
      real(kind=ISFFT_REAL_KIND), allocatable :: r(:, :, :)

      call Park_Field(field)
      call Get_Real(field, shape(input), r)
      call move_alloc(r, field%r)
      field%r = input
   end subroutine ! Load_Field

   subroutine Destroy_Poisson_Operator(op)
      type(Poisson_Operator), intent(inout) :: op

      if (allocated(op%multiplier)) deallocate (op%multiplier)
   end subroutine ! Destroy_Poisson_Operator

   !------------------------------------------------------------------------------
   !! Shared transport workspace lifecycle
   !------------------------------------------------------------------------------
   !! Seal only after every selected execution path has reserved its scratch.
   subroutine Seal_Layout_Workspace()

      workspace_sealed = .true.
   end subroutine ! Seal_Layout_Workspace

   subroutine Reserve_R(storage, required_capacity)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: storage(:)
      integer(ISFFT_INT64_KIND), intent(in) :: required_capacity

      if (allocated(storage)) then
         if (size(storage, kind=ISFFT_INT64_KIND) >= required_capacity) return
         if (workspace_sealed) error stop 'real MPI workspace was not reserved'
         deallocate (storage)
      end if
      if (workspace_sealed) error stop 'real MPI workspace was not reserved'
      allocate (storage(required_capacity))
   end subroutine ! Reserve_R

   subroutine Reserve_Z(storage, required_capacity)
      complex(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: storage(:)
      integer(ISFFT_INT64_KIND), intent(in) :: required_capacity

      if (allocated(storage)) then
         if (size(storage, kind=ISFFT_INT64_KIND) >= required_capacity) return
         if (workspace_sealed) error stop 'complex MPI workspace was not reserved'
         deallocate (storage)
      end if
      if (workspace_sealed) error stop 'complex MPI workspace was not reserved'
      allocate (storage(required_capacity))
   end subroutine ! Reserve_Z

   subroutine Release_Layout_Workspace()

      workspace_sealed = .false.
      if (allocated(send_real_buffer)) deallocate (send_real_buffer)
      if (allocated(recv_real_buffer)) deallocate (recv_real_buffer)
      if (allocated(send_complex_buffer)) deallocate (send_complex_buffer)
      if (allocated(recv_complex_buffer)) deallocate (recv_complex_buffer)
   end subroutine ! Release_Layout_Workspace

   !------------------------------------------------------------------------------
   !! Bounded-memory streaming and chunk packing
   !------------------------------------------------------------------------------
   !! Keep large-count scratch storage bounded even after earlier ordinary transfers.
   subroutine Bound_Layout_Workspace(capacity)
      integer(ISFFT_INT64_KIND), intent(in) :: capacity

      if (workspace_sealed) return
      if (allocated(send_real_buffer)) then
         if (size(send_real_buffer, kind=ISFFT_INT64_KIND) > capacity) deallocate (send_real_buffer)
      end if
      if (allocated(recv_real_buffer)) then
         if (size(recv_real_buffer, kind=ISFFT_INT64_KIND) > capacity) deallocate (recv_real_buffer)
      end if
      if (allocated(send_complex_buffer)) then
         if (size(send_complex_buffer, kind=ISFFT_INT64_KIND) > capacity) deallocate (send_complex_buffer)
      end if
      if (allocated(recv_complex_buffer)) then
         if (size(recv_complex_buffer, kind=ISFFT_INT64_KIND) > capacity) deallocate (recv_complex_buffer)
      end if
   end subroutine ! Bound_Layout_Workspace

end module ISFFT_Workspace
