!------------------------------------------------------------------------------
! ISFFT_Planning
! Estimate rules and Measure selection using the shared execution paths.
!------------------------------------------------------------------------------
module ISFFT_Planning
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND, ISFFT_REAL_MPI_TYPE
   use ISFFT_Plans, only: Axis_Plan, axis_plan_cache, block_cache, Block_Choice, &
      Kernel_Plan, poisson_plan, prepared_kernel_plans, public_export, &
      public_import, public_plans, redistribution_forward_plan, redistribution_inverse_plan, &
      Transfer_Plan, Transform_Plan
   use ISFFT_Problem_State, only: internal_layout, ISFFT_PLAN_ESTIMATE, ISFFT_PLAN_MEASURE, Layout, &
      layout_policy, original_layout, planning_level, prepare_scope, &
      public_fft_modes, public_layout, topology_ready, transport_policy
   use ISFFT_Workspace, only: Access_Candidate, Clear_Field, fft_precompute_slots, poisson_symbol, &
      Release_Layout_Workspace, Require_FFT_Precompute, sample_cache, Sample_Scores, &
      Spectral_Field
   use ISFFT_MPI_Part, only: Axis_Position, Layout_Box, Layout_Shape, Make_Layout, Same_Layout
   use ISFFT_DCT_kernel, only: ISFFT_Execute_DCT_Plan
   use ISFFT_DST_kernel, only: ISFFT_Execute_DST_Plan
   use ISFFT_FFT_kernel, only: ISFFT_Execute_FFT_Complex_Strategy, ISFFT_Execute_R2C_Batch, &
      ISFFT_Execute_R2C_Generic_Batch, ISFFT_Execute_R2C_Odd_Batch, &
      ISFFT_Execute_R2C_Scalar, ISFFT_Is_Three_Times_Power_Of_Two
   use ISFFT_Axis_Execution, only: Execute_Kernel
   use ISFFT_Communication, only: Build_Block_Gather, Build_Transfer, datatype_a2a, datatype_p2p, &
      memorysave_p2p, Move_Field, packed_a2a, packed_p2p, &
      Release_Transform_Plan, resources, Select_Transfer_Method
   use ISFFT_Transform_Execution, only: Block_Sample_Lines, Execute_Block_Axis
   use ISFFT_Allocation, only: Finalize_Plans, Reserve_Axis_Scratch, Reserve_Block, Reserve_FFT_Axis, &
      Reserve_FFT_Base, Reserve_Field, Reserve_R2C_Workspace
   use ISFFT_Precompute, only: Build_Poisson_Operator, Ensure_FFT_Precompute, Ensure_R2C_Precompute, &
      Ensure_Real_Transform_Precompute, &
      Prepare_Complex_Candidates, Prepare_R2C_Candidates, Prepare_Real_Workspace
   implicit none
   private
   public :: Build_Kernel_Plan, Build_Transform_Plan, COMMUNICATION_TRIAL, Ensure_Public_Plan, Estimate_Chunk
   public :: Estimate_Width, FULL_BATCH_TRIAL, Initialize_Plans, ISFFT_Clone_Complex_Strategy, ISFFT_Clone_R2C_Strategy
   public :: ISFFT_Forget_Complex_Strategy, ISFFT_Forget_R2C_Strategy, ISFFT_Prepare_Kernel, &
      ISFFT_Store_Complex_Strategy, ISFFT_Store_R2C_Strategy
   public :: ISFFT_Tune_FFT_Complex_Batch, ISFFT_Tune_FFT_Complex_Shape, &
      ISFFT_Tune_FFT_R2C_Batch, ISFFT_Tune_FFT_R2C_Shape, KERNEL_TRIAL
   public :: Plan_Communication_Strategy, Plan_Transport, Policy_Axis_Start, Policy_Axis_Widths, Policy_Block_Capacity
   public :: Policy_Block_Start, Policy_Configure, Policy_Repeat_Count, Policy_Repeats, Policy_Sample
   public :: Policy_Warmup_Count, Policy_Warmups, Policy_Widths, SAMPLED_BATCH_TRIAL, Select_Best_Candidate
   public :: Tune_Access_And_Batching, Tune_Block_Width, Tune_Real_Kernel_And_Inner_Batch, Tune_Transfer

   !! Initialization timing stages.
   integer, parameter :: KERNEL_TRIAL = 1
   integer, parameter :: SAMPLED_BATCH_TRIAL = 2
   integer, parameter :: COMMUNICATION_TRIAL = 3
   integer, parameter :: FULL_BATCH_TRIAL = 4

   !! Working-set budgets in bytes.
   integer(ISFFT_INT64_KIND), parameter :: default_cache_budget_bytes = 131072_ISFFT_INT64_KIND !! 128 KiB
   integer(ISFFT_INT64_KIND), parameter :: strided_cache_budget_bytes = 262144_ISFFT_INT64_KIND !! 256 KiB
   integer(ISFFT_INT64_KIND), parameter :: min_cache_budget_bytes = 16384_ISFFT_INT64_KIND !! 16 KiB
   integer(ISFFT_INT64_KIND), parameter :: max_cache_budget_bytes = 1048576_ISFFT_INT64_KIND !! 1 MiB
   integer(ISFFT_INT64_KIND), parameter :: stream_workspace_budget_bytes = 8388608_ISFFT_INT64_KIND !! 8 MiB

   integer, parameter :: Policy_Warmups = 2
   integer, save :: node_ranks = 1
   integer(ISFFT_INT64_KIND), save :: cache_budget = 0
   !! Keep level-indexed tables aligned with ESTIMATE=0 and MEASURE=1.
   integer, parameter :: ordinary_repeats(0:1) = [0, 3], communication_repeats(0:1) = [0, 2]
   integer(ISFFT_INT64_KIND), parameter :: max_packed_workspace_bytes = 2147483648_ISFFT_INT64_KIND

   integer, save :: plan_dim, plan_bc(3), plan_decomp, plan_generation

contains

   !------------------------------------------------------------------------------
   !! Selection policy, candidate widths and measurement counts
   !------------------------------------------------------------------------------
   integer function Policy_Warmup_Count(plan_level, trial_stage) result(n)
      integer, intent(in) :: plan_level, trial_stage

      n = 2
      if (plan_level == 0) n = 0
      if (plan_level == 1 .and. trial_stage /= KERNEL_TRIAL) n = 1
   end function ! Policy_Warmup_Count

   integer function Policy_Repeat_Count(plan_level, trial_stage) result(n)
      integer, intent(in) :: plan_level, trial_stage

      n = Policy_Repeats(plan_level, trial_stage == COMMUNICATION_TRIAL)
      if (plan_level == 1 .and. trial_stage == SAMPLED_BATCH_TRIAL) n = 2
      if (plan_level == 1 .and. trial_stage == FULL_BATCH_TRIAL) n = 1
   end function ! Policy_Repeat_Count

   subroutine Policy_Configure(plan_level, collect_environment)
      integer, intent(in) :: plan_level
      logical, intent(in), optional :: collect_environment
      logical :: use_environment
      integer :: i, ios, rank, ierr, shared, u, cpu, idx, lev, socket_id, core_id
      integer(ISFFT_INT64_KIND) :: bytes, minbudget
      character(2048) :: path, line, allowed, shared_list
      character(64) :: v
      character(24) :: cs, ix

      if (plan_level < 0 .or. plan_level > 1) error stop "plan_level must be ESTIMATE (0) or MEASURE (1)"
      use_environment = .true.
      if (present(collect_environment)) use_environment = collect_environment
      cache_budget = 0
      node_ranks = 1
      allowed = 'unknown'
      shared_list = 'unknown'
      cpu = 0
      socket_id = -1
      core_id = -1
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      if (use_environment) then
         call MPI_Comm_split_type(ISFFT_Comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, shared, ierr)
         if (ierr == MPI_SUCCESS) then
            call MPI_Comm_size(shared, node_ranks, ierr)
            call MPI_Comm_free(shared, ierr)
         end if
         !! Launcher can supply read-only macOS sysctl data. No external command is
         !! executed by the library; absent information always falls back to a rule.
         call get_environment_variable('ISFFT_CACHE_BYTES', v, status=ios)
         if (ios == 0) then
            bytes = 0
            read (v, *, iostat=ios) bytes
            if (ios == 0 .and. bytes > 0) cache_budget = bytes / max(1, node_ranks) / 4
         end if
         open (newunit=u, file='/proc/self/status', status='old', iostat=ios)
         if (ios == 0) then
            do
               read (u, '(A)', iostat=ios) line
               if (ios /= 0) exit
               if (index(line, 'Cpus_allowed_list:') == 1) then
                  allowed = adjustl(line(19:))
                  v = allowed
                  do i = 1, len_trim(v)
                     if (v(i:i) == '-' .or. v(i:i) == ',') v(i:i) = ' '
                  end do
                  read (v, *, iostat=ios) cpu
                  exit
               end if
            end do
            close (u)
         end if
         write (cs, '(I0)') cpu
         path = '/sys/devices/system/cpu/cpu'//trim(cs)//'/topology/'
         open (newunit=u, file=trim(path)//'physical_package_id', status='old', iostat=ios)
         if (ios == 0) then
            read (u, *, iostat=ios) socket_id
            close (u)
         end if
         open (newunit=u, file=trim(path)//'core_id', status='old', iostat=ios)
         if (ios == 0) then
            read (u, *, iostat=ios) core_id
            close (u)
         end if
         do idx = 0, 7
            write (ix, '(I0)') idx
            path = '/sys/devices/system/cpu/cpu'//trim(cs)//'/cache/index'//trim(ix)//'/'
            open (newunit=u, file=trim(path)//'level', status='old', iostat=ios)
            if (ios /= 0) cycle
            read (u, *, iostat=ios) lev
            close (u)
            if (ios /= 0 .or. lev < 2) cycle
            open (newunit=u, file=trim(path)//'size', status='old', iostat=ios)
            if (ios /= 0) cycle
            read (u, '(A)', iostat=ios) v
            close (u)
            if (ios /= 0) cycle
            i = len_trim(v)
            bytes = 0
            if (i > 1) then
               read (v(:i - 1), *, iostat=ios) bytes
               if (ios /= 0) cycle
               if (v(i:i) == 'K') bytes = bytes * 1024
               if (v(i:i) == 'M') bytes = bytes * 1024 * 1024
            end if
            open (newunit=u, file=trim(path)//'shared_cpu_list', status='old', iostat=ios)
            if (ios == 0) then
               read (u, '(A)', iostat=ios) shared_list
               close (u)
            end if
            !! Conservative share: divide by all ranks on this node. Never treat a
            !! shared cache as rank-private; affinity migration makes exact shares unknown.
            cache_budget = max(cache_budget, bytes / max(1, node_ranks) / 4)
         end do
         !! A single candidate order must be identical on every rank.
         call MPI_Allreduce(cache_budget, minbudget, 1, MPI_INTEGER8, MPI_MIN, ISFFT_Comm, ierr)
         cache_budget = minbudget
      end if
   end subroutine ! Policy_Configure

   subroutine Policy_Axis_Start(plan_level, transform_length, transform_kind, local_shape, memory_axis, &
      complex_input, outer_batch_width, direct_access, direct_batch_width)
      integer, intent(in) :: plan_level, transform_length, transform_kind, local_shape(3), memory_axis
      logical, intent(in) :: complex_input
      integer, intent(out) :: outer_batch_width, direct_batch_width
      logical, intent(out) :: direct_access
      integer :: lines, bytes
      integer(ISFFT_INT64_KIND) :: budget, footprint, stride, target

      lines = max(1, int(product(int(local_shape, ISFFT_INT64_KIND)) / max(1, transform_length)))
      bytes = merge(16, 8, transform_kind == 2 .or. complex_input)
      outer_batch_width = Estimate_Width(transform_length, lines, bytes, 3, lines, memory_axis)
      direct_access = .false.
      direct_batch_width = 1
      if (memory_axis == 1) outer_batch_width = lines
      stride = 1
      if (memory_axis > 1) stride = product(int(local_shape(:memory_axis - 1), ISFFT_INT64_KIND))
      budget = cache_budget
      if (budget <= 0) budget = merge(strided_cache_budget_bytes, default_cache_budget_bytes, stride > 1)
      budget = max(min_cache_budget_bytes, min(max_cache_budget_bytes, budget))
      footprint = max(1_ISFFT_INT64_KIND, int(transform_length, ISFFT_INT64_KIND) * bytes * 3)
      target = max(1_ISFFT_INT64_KIND, min(int(lines, ISFFT_INT64_KIND), budget / footprint))
      outer_batch_width = 1
      do while (outer_batch_width <= target / 2)
         outer_batch_width = outer_batch_width * 2
      end do
      if (memory_axis == 1 .and. (transform_kind == 1 .or. transform_kind >= 3)) outer_batch_width = lines
      !! The existing direct executor supports all prepared C2C algorithms and tails.
      if (transform_kind == 2 .and. memory_axis > 1) then
         direct_access = .true.
         do while (direct_batch_width <= min(16, lines) / 2)
            direct_batch_width = direct_batch_width * 2
         end do
         outer_batch_width = 1
      end if
   end subroutine ! Policy_Axis_Start

   subroutine Policy_Axis_Widths(plan_level, estimated_width, max_batch_width, candidate_widths, transform_length, &
      bytes_per_element, memory_axis)
      integer, intent(in) :: plan_level, estimated_width, max_batch_width, transform_length, bytes_per_element, &
         memory_axis
      integer, allocatable, intent(out) :: candidate_widths(:)
      integer :: a

      call Policy_Widths(plan_level, estimated_width, max_batch_width, 4, candidate_widths)
      if (plan_level /= 1) return
      a = Estimate_Width(transform_length, max_batch_width, bytes_per_element, 3, max_batch_width, memory_axis)
   end subroutine ! Policy_Axis_Widths

   integer function Policy_Block_Capacity(transform_length, line_count) result(width)
      integer, intent(in) :: transform_length, line_count
      integer(ISFFT_INT64_KIND), parameter :: block_buffer_budget_bytes = 1048576_ISFFT_INT64_KIND !! 1 MiB

      width = int(max(1_ISFFT_INT64_KIND, &
         min(int(max(1, line_count), ISFFT_INT64_KIND), block_buffer_budget_bytes / &
            (16_ISFFT_INT64_KIND * max(1, transform_length)))))
   end function ! Policy_Block_Capacity

   integer function Policy_Block_Start(transform_length, line_count, bytes_per_element, axis_process_count, &
      batch_capacity) result(width)
      integer, intent(in) :: transform_length, line_count, bytes_per_element, axis_process_count, batch_capacity
      integer(ISFFT_INT64_KIND) :: budget, footprint, target

      !! Gathered line, packed messages and transform scratch share the cache.
      !! Keep the execution tile capacity, and account for local message storage.
      budget = cache_budget
      if (budget <= 0) budget = default_cache_budget_bytes
      budget = max(min_cache_budget_bytes, min(max_cache_budget_bytes, budget))
      footprint = int(bytes_per_element, ISFFT_INT64_KIND) * &
         (3_ISFFT_INT64_KIND * max(1, transform_length) + &
         (int(max(1, transform_length), ISFFT_INT64_KIND) + max(1, axis_process_count) - 1) / &
            max(1, axis_process_count))
      target = max(1_ISFFT_INT64_KIND, &
         min(int(max(1, min(line_count, batch_capacity)), ISFFT_INT64_KIND), budget / &
            max(1_ISFFT_INT64_KIND, footprint)))
      width = 1
      do while (width <= target / 2)
         width = 2 * width
      end do
   end function ! Policy_Block_Start

   integer function Select_Best_Candidate(candidate_times, baseline_index, timer_resolution) result(best)
      real(ISFFT_REAL_KIND), intent(in) :: candidate_times(:), timer_resolution
      integer, intent(in) :: baseline_index

      best = minloc(candidate_times, dim=1)
      !! A benefit smaller than one timer tick is not resolvable. Retain the
      !! estimated legal candidate; never extend the repetition count implicitly.
      if (candidate_times(baseline_index) - candidate_times(best) <= timer_resolution) then
         best = baseline_index
      end if
   end function ! Select_Best_Candidate

   pure integer function Policy_Repeats(plan_level, communication) result(n)
      integer, intent(in) :: plan_level
      logical, intent(in), optional :: communication

      n = ordinary_repeats(max(0, min(1, plan_level)))
      if (present(communication)) then
         if (communication) n = communication_repeats(max(0, min(1, plan_level)))
      end if
   end function ! Policy_Repeats

   pure integer function Policy_Sample(plan_level, line_count, access_trial) result(n)
      integer, intent(in) :: plan_level, line_count
      logical, intent(in), optional :: access_trial
      integer :: limit

      limit = 64
      if (present(access_trial)) then
         if (access_trial) limit = 4 * limit
      end if
      n = max(1, min(line_count, limit))
   end function ! Policy_Sample

   pure integer function Estimate_Width(transform_length, line_count, bytes_per_element, workspace_copies, &
      batch_capacity, memory_axis) result(w)
      integer, intent(in) :: transform_length, line_count, bytes_per_element, workspace_copies, batch_capacity
      integer, intent(in), optional :: memory_axis
      integer(ISFFT_INT64_KIND) :: footprint, budget
      integer :: target

      !! 128 KiB is a portable working-set starting point, not a detected cache size.
      !! Packing/strided work includes both gathered values and the transform workspace.
      budget = default_cache_budget_bytes
      if (present(memory_axis)) then
         if (memory_axis > 1) budget = strided_cache_budget_bytes
      end if
      footprint = max(1_ISFFT_INT64_KIND, int(transform_length, ISFFT_INT64_KIND) * bytes_per_element * &
         workspace_copies)
      target = int(max(1_ISFFT_INT64_KIND, &
         min(int(min(line_count, batch_capacity), ISFFT_INT64_KIND), budget / footprint)))
      w = 1
      do while (w <= target / 2)
         w = 2 * w
      end do
   end function ! Estimate_Width

   pure integer(ISFFT_INT64_KIND) function Estimate_Chunk(elements, bytes_per_element) result(w)
      integer(ISFFT_INT64_KIND), intent(in) :: elements
      integer, intent(in) :: bytes_per_element
      integer(ISFFT_INT64_KIND) :: target

      !! Two streaming buffers share an 8 MiB working-set budget.
      target = max(1_ISFFT_INT64_KIND, min(elements, stream_workspace_budget_bytes / max(1, 2 * bytes_per_element)))
      w = 1
      do while (w <= target / 2)
         w = 2 * w
      end do
   end function ! Estimate_Chunk

   subroutine Policy_Widths(plan_level, estimated_width, max_batch_width, candidate_family, candidate_widths)
      integer, intent(in) :: plan_level, estimated_width, max_batch_width, candidate_family
      integer, allocatable, intent(out) :: candidate_widths(:)
      integer :: values(512), count, i, first
      integer :: base(8)
      integer, parameter :: fft_widths(8) = [1, 2, 3, 4, 6, 8, 12, 16]
      integer, parameter :: real_widths(8) = [1, 4, 8, 16, 24, 32, 40, 64]
      !! Odd R2C widths count paired lines, not internal FFT batch capacity.
      integer, parameter :: paired_widths(8) = [1, 2, 4, 8, 16, 24, 32, 64]

      count = 0
      if (candidate_family <= 3) then
         select case (candidate_family)
         case (1)
            base = fft_widths
         case (2)
            base = real_widths
         case (3)
            base = paired_widths
         end select
         !! Start at the nearest legal template point; do not silently add a ninth
         !! MEASURE width. ESTIMATE itself retains its unchanged empirical rule.
         first = 0
         do i = 1, size(base)
            if (base(i) > max_batch_width) cycle
            if (first == 0) then
               first = i
            else if (abs(base(i) - estimated_width) < abs(base(first) - estimated_width)) then
               first = i
            end if
         end do
         if (first > 0) call add(base(first))
         do i = 1, size(base)
            call add(base(i))
         end do
      else
         call add(estimated_width)
         call add(max(1, estimated_width / 2))
         call add(min(max_batch_width, 2 * estimated_width))
         do i = 0, 8
            call add(2**i)
         end do
         call add(max_batch_width)
      end if
      candidate_widths = values(:count)
   contains

      subroutine add(value)
         integer, intent(in) :: value

         if (value < 1 .or. value > max_batch_width) return
         if (any(values(:count) == value)) return
         count = count + 1
         if (count > size(values)) error stop 'policy candidate capacity exceeded'
         values(count) = value
      end subroutine ! add
   end subroutine ! Policy_Widths

   !------------------------------------------------------------------------------
   !! Communication strategy selection
   !------------------------------------------------------------------------------
   subroutine Tune_Transfer(plan, level)
      type(Transfer_Plan), intent(inout) :: plan
      integer, intent(in), optional :: level
      integer :: effort, id

      effort = 1
      if (present(level)) effort = level
      if (effort < 0 .or. effort > 1) error stop 'Invalid transport planning level'
      if (plan%local_only) return
      id = plan%resource_id
      if (id > 0) then
         if (resources(id)%tuned(effort)) then
            plan%direct = resources(id)%direct(effort)
            plan%streamed = resources(id)%streamed(effort)
            plan%packed_p2p_active = resources(id)%packed_p2p_active(effort)
            plan%stream_width = resources(id)%width(effort)
            plan%sparse = resources(id)%sparse(effort)
            return
         end if
      end if
      call Select_Communication_Strategy(plan, effort)
      if (id > 0) then
         resources(id)%tuned(effort) = .true.
         resources(id)%direct(effort) = plan%direct
         resources(id)%streamed(effort) = plan%streamed
         resources(id)%width(effort) = plan%stream_width
         resources(id)%sparse(effort) = plan%sparse
         resources(id)%packed_p2p_active(effort) = plan%packed_p2p_active
      end if
   end subroutine ! Tune_Transfer

   !! All ranks must visit candidates and collectives in the same order.
   !! Keep coordinate checks coupled to each measured candidate.
   !! Tune transport independently of FFT kernel widths and axis batch widths.
   !! level: 0=ESTIMATE, 1=MEASURE. Initialization is collective.
   subroutine Select_Communication_Strategy(plan, level)
      type(Transfer_Plan), intent(inout) :: plan
      integer, intent(in), optional :: level
      real(kind=ISFFT_REAL_KIND), allocatable :: real_input(:, :, :), real_output(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: complex_input(:, :, :), complex_output(:, :, :)
      integer :: rank, ierr, local_shape(3), candidate_index, repeat_index, candidate_count, repeat_count, effort, &
         best_candidate
      integer :: candidate_methods(24), local_bad, global_bad, alloc_status
      integer :: warmup_count, timed_count, check_count
      integer(ISFFT_INT64_KIND) :: candidate_widths(24), stream_width_choices(20), candidate_width, &
         candidate_workspace_elements(24), local_workspace_elements, packed_workspace_elements, bytes_per_element
      real(kind=ISFFT_REAL_KIND) :: start_time, elapsed, global_elapsed, score, best_score, scores(24)
      logical :: light

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      effort = 1
      if (present(level)) effort = level
      if (effort < 0 .or. effort > 1) error stop 'Invalid transport planning level'
      plan%packed_p2p_active = .false.
      plan%sparse = .false.
      plan%direct = .false.
      plan%streamed = plan%chunked
      plan%stream_width = &
         Estimate_Chunk(max(1_ISFFT_INT64_KIND, min(plan%chunk_limit, plan%max_peer)), &
            merge(16, 8, plan%source%complex_data))
      if (effort == 0) then
         plan%direct = .not. plan%chunked .and. plan%gather_comm == MPI_COMM_NULL
         return
      end if
      if (plan%local_only) return
      warmup_count = 0
      timed_count = 0
      check_count = 0
      light = .false.
      candidate_count = 0
      candidate_methods = 0
      candidate_widths = 0
      candidate_workspace_elements = 0
      if (.not. plan%chunked) then
         !! MPI datatypes avoid the library's pack buffers for non-replicated exchanges.
         if (plan%gather_comm == MPI_COMM_NULL) then
            candidate_count = candidate_count + 1
            candidate_methods(candidate_count) = datatype_a2a
            candidate_count = candidate_count + 1
            candidate_methods(candidate_count) = datatype_p2p
         end if
         local_workspace_elements = sum(plan%send_counts) + sum(plan%recv_counts)
         if (plan%gather_comm /= MPI_COMM_NULL) local_workspace_elements = plan%gather_count + sum(plan%recv_counts)
         call MPI_Allreduce(local_workspace_elements, packed_workspace_elements, 1, MPI_INTEGER8, MPI_MAX, &
            ISFFT_Comm, ierr)
         bytes_per_element = storage_size(0.0_ISFFT_REAL_KIND) / 8
         if (plan%source%complex_data) bytes_per_element = 2 * bytes_per_element
         !! Bound the optional full-pack candidate to 2 GiB of extra library scratch per rank.
         !! This does not bound user fields or internal storage allocated by MPI itself.
         if (packed_workspace_elements <= max_packed_workspace_bytes / bytes_per_element) then
            candidate_count = candidate_count + 1
            candidate_methods(candidate_count) = packed_a2a
            candidate_workspace_elements(candidate_count) = packed_workspace_elements
            if (plan%gather_comm == MPI_COMM_NULL .and. allocated(plan%requests)) then
               candidate_count = candidate_count + 1
               candidate_methods(candidate_count) = packed_p2p
               candidate_workspace_elements(candidate_count) = packed_workspace_elements
            end if
         end if
      end if
      stream_width_choices = 0
      stream_width_choices(1:3) = &
         [plan%stream_width, max(1_ISFFT_INT64_KIND, plan%stream_width / 2), 2 * plan%stream_width]
      repeat_count = Policy_Repeats(effort, .true.)
      do candidate_index = 1, size(stream_width_choices)
         if (stream_width_choices(candidate_index) == 0) cycle
         candidate_width = &
            max(1_ISFFT_INT64_KIND, min(stream_width_choices(candidate_index), plan%chunk_limit, plan%max_peer))
         if (any(candidate_widths(:candidate_count) == candidate_width)) cycle
         candidate_count = candidate_count + 1
         candidate_methods(candidate_count) = memorysave_p2p
         candidate_widths(candidate_count) = candidate_width
         candidate_workspace_elements(candidate_count) = 2 * candidate_width
      end do
      !! Only one input and one output are needed. Check the known global-coordinate
      !! pattern directly, rather than retaining an additional full reference field.
      local_shape = Layout_Shape(plan%source, rank)
      if (plan%source%complex_data) then
         allocate (complex_input(local_shape(1), local_shape(2), local_shape(3)), stat=alloc_status)
      else
         allocate (real_input(local_shape(1), local_shape(2), local_shape(3)), stat=alloc_status)
      end if
      local_bad = merge(1, 0, alloc_status /= 0)
      call MPI_Allreduce(local_bad, global_bad, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      if (global_bad /= 0) error stop 'Cannot allocate transport planning input'
      call Fill_Sample()
      local_shape = Layout_Shape(plan%destination, rank)
      if (plan%source%complex_data) then
         allocate (complex_output(local_shape(1), local_shape(2), local_shape(3)), stat=alloc_status)
      else
         allocate (real_output(local_shape(1), local_shape(2), local_shape(3)), stat=alloc_status)
      end if
      local_bad = merge(1, 0, alloc_status /= 0)
      call MPI_Allreduce(local_bad, global_bad, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      if (global_bad /= 0) error stop 'Cannot allocate transport planning output'
      best_score = huge(best_score)
      best_candidate = 1
      scores = huge(best_score)
      do candidate_index = 1, candidate_count
         call Select_Candidate(candidate_index)
         !! Remove previous candidate buffers before warming the new path.
         call Release_Layout_Workspace()
         !! A previous candidate must not hide a missing write in this candidate.
         if (plan%source%complex_data) then
            complex_output = cmplx(-huge(0.0_ISFFT_REAL_KIND), 0.0_ISFFT_REAL_KIND, ISFFT_REAL_KIND)
         else
            real_output = -huge(0.0_ISFFT_REAL_KIND)
         end if
         do repeat_index = 1, Policy_Warmup_Count(effort, COMMUNICATION_TRIAL)
            start_time = MPI_Wtime()
            call Execute_Sample()
            warmup_count = warmup_count + 1
            elapsed = MPI_Wtime() - start_time
         end do
         score = 0
         do repeat_index = 1, repeat_count
            call MPI_Barrier(ISFFT_Comm, ierr)
            start_time = MPI_Wtime()
            call Execute_Sample()
            elapsed = MPI_Wtime() - start_time
            timed_count = timed_count + 1
            call MPI_Allreduce(elapsed, global_elapsed, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
            score = score + global_elapsed
         end do
         call Check_Sample()
         score = score / repeat_count
         scores(candidate_index) = score

         !! Near ties favor less library scratch storage. This also avoids promoting
         !! a noisily faster streaming candidate over a nearly identical direct path.
         if (candidate_index == 1) then
            best_score = score
            best_candidate = candidate_index
         else if (score < best_score) then
            best_score = score
            best_candidate = candidate_index
         else if (score == best_score .and. &
            candidate_workspace_elements(candidate_index) < candidate_workspace_elements(best_candidate)) then
            best_score = score
            best_candidate = candidate_index
         end if
      end do
      call Select_Candidate(best_candidate)
      call Release_Layout_Workspace()
   contains

      subroutine Select_Candidate(candidate_index)
         integer, intent(in) :: candidate_index

         call Select_Transfer_Method(plan, candidate_methods(candidate_index))
         if (plan%streamed) then
            plan%stream_width = candidate_widths(candidate_index)
         else
            plan%stream_width = Estimate_Chunk( &
               max(1_ISFFT_INT64_KIND, min(plan%chunk_limit, plan%max_peer)), merge(16, 8, plan%source%complex_data))
         end if
      end subroutine ! Select_Candidate

      subroutine Execute_Sample()

         if (plan%source%complex_data) then
            call Move_Field(plan, complex_input, complex_output)
         else
            call Move_Field(plan, real_input, real_output)
         end if
      end subroutine ! Execute_Sample

      subroutine Check_Sample()
         integer :: i, j, k, physical(3), offset(3), dims(3)
         integer(ISFFT_INT64_KIND) :: stride(3), step(3), base
         real(kind=ISFFT_REAL_KIND) :: value

         call Layout_Box(plan%destination, rank, physical, offset)
         dims = physical(plan%destination%axis_order)
         stride = [1_ISFFT_INT64_KIND, int(plan%destination%global_shape(1), ISFFT_INT64_KIND), &
            int(plan%destination%global_shape(1), ISFFT_INT64_KIND) * plan%destination%global_shape(2)]
         step = stride(plan%destination%axis_order)
         local_bad = 0
         do k = 1, dims(3)
            do j = 1, dims(2)
               base = sum(int(offset, ISFFT_INT64_KIND) * stride) + int(j - 1, ISFFT_INT64_KIND) * step(2) + &
                  int(k - 1, ISFFT_INT64_KIND) * step(3)
               do i = 1, dims(1)
                  value = real(base + int(i - 1, ISFFT_INT64_KIND) * step(1), ISFFT_REAL_KIND)
                  if (plan%source%complex_data) then
                     if (complex_output(i, j, k) /= cmplx(value, -0.5_ISFFT_REAL_KIND * value, ISFFT_REAL_KIND)) &
                        local_bad = 1
                  else
                     if (real_output(i, j, k) /= value) local_bad = 1
                  end if
               end do
            end do
         end do
         call MPI_Allreduce(local_bad, global_bad, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
         if (global_bad /= 0) error stop 'transport candidate failed coordinate check'
      end subroutine ! Check_Sample

      subroutine Fill_Sample()
         integer :: i, j, k, physical(3), offset(3), dims(3)
         integer(ISFFT_INT64_KIND) :: stride(3), step(3), base
         real(kind=ISFFT_REAL_KIND) :: value

         call Layout_Box(plan%source, rank, physical, offset)
         dims = physical(plan%source%axis_order)
         stride = [1_ISFFT_INT64_KIND, int(plan%source%global_shape(1), ISFFT_INT64_KIND), &
            int(plan%source%global_shape(1), ISFFT_INT64_KIND) * &
            plan%source%global_shape(2)]
         step = stride(plan%source%axis_order)
         do k = 1, dims(3)
            do j = 1, dims(2)
               base = sum(int(offset, ISFFT_INT64_KIND) * stride) + int(j - 1, ISFFT_INT64_KIND) * step(2) + &
                  int(k - 1, ISFFT_INT64_KIND) * step(3)
               do i = 1, dims(1)
                  value = real(base + int(i - 1, ISFFT_INT64_KIND) * step(1), ISFFT_REAL_KIND)
                  if (plan%source%complex_data) then
                     complex_input(i, j, k) = cmplx(value, -0.5_ISFFT_REAL_KIND * value, ISFFT_REAL_KIND)
                  else
                     real_input(i, j, k) = value
                  end if
               end do
            end do
         end do
      end subroutine ! Fill_Sample

   end subroutine ! Select_Communication_Strategy

   !------------------------------------------------------------------------------
   !! Numerical kernel and inner-batch selection
   !------------------------------------------------------------------------------

   !! Low-level initialization for direct kernel users. n is the full grid length:
   !! kind 1=R2C, 2=C2C, 3=DST-I (n-2 interior values), 4=DCT-I (n values).
   !! MPI must be initialized when tune is true. High-level users use ISFFT_Initialize.
   recursive subroutine ISFFT_Prepare_Kernel(transform_length, transform_kind, line_count, tune, plan, axis_storage)
      integer, intent(in) :: transform_length, transform_kind, line_count
      logical, intent(in), optional :: tune, axis_storage
      type(Kernel_Plan), intent(out), optional :: plan

      call ISFFT_Prepare_Kernel_Impl(transform_length, transform_kind, line_count, tune, plan, axis_storage)
   end subroutine ! ISFFT_Prepare_Kernel

   !! A local cache hit still participates when another rank needs collective tuning.
   recursive subroutine ISFFT_Prepare_Kernel_Impl(transform_length, transform_kind, line_count, tune, plan, &
      axis_storage)
      integer, intent(in) :: transform_length, transform_kind, line_count
      logical, intent(in), optional :: tune, axis_storage
      type(Kernel_Plan), intent(out), optional :: plan
      integer :: slot, l, e, count, hit, all_hit, ierr
      type(Kernel_Plan), allocatable :: grown(:)
      logical :: do_tune
      type(Kernel_Plan) :: p

      if (transform_length < 1 .or. line_count < 1 .or. transform_kind < 1 .or. transform_kind > 4) &
         error stop 'ISFFT_Prepare_Kernel: invalid shape/kind'
      do_tune = .true.
      if (present(tune)) do_tune = tune
      hit = 0
      if (do_tune .and. planning_level /= ISFFT_PLAN_ESTIMATE) then
         hit = 0
         if (allocated(prepared_kernel_plans)) then
            do e = 1, size(prepared_kernel_plans)
               if (prepared_kernel_plans(e)%transform_length == transform_length .and. &
                  prepared_kernel_plans(e)%transform_kind == transform_kind .and. &
                  prepared_kernel_plans(e)%line_count == line_count) then
                  hit = e
                  exit
               end if
            end do
         end if
         e = 0
         if (hit > 0) e = 1
         call MPI_Allreduce(e, all_hit, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         if (all_hit == 1) then
            p = prepared_kernel_plans(hit)
            if (transform_kind == 2) then
               if (present(axis_storage)) then
                  if (axis_storage) call Reserve_FFT_Axis(transform_length)
               else
                  call Reserve_FFT_Axis(transform_length)
               end if
            end if
            if (present(plan)) plan = p
            return
         end if
      end if
      if (transform_kind <= 2) then
         if (transform_kind == 1) then
            call Ensure_R2C_Precompute(transform_length, slot)
            call Reserve_R2C_Workspace(transform_length, line_count)
         else
            call Ensure_FFT_Precompute(transform_length, slot)
         end if
      else
         l = max(1, transform_length - 1)
         call Ensure_Real_Transform_Precompute(l, slot)
         call Prepare_Real_Workspace(l, line_count)
      end if
      if (transform_kind == 2) then
         if (present(axis_storage)) then
            if (axis_storage) then
               call Reserve_FFT_Axis(transform_length)
            end if
         else
            call Reserve_FFT_Axis(transform_length)
         end if
      end if
      p%transform_length = transform_length
      p%transform_kind = transform_kind
      p%line_count = line_count
      p%outer_batch_width = line_count
      if (.not. do_tune .or. planning_level == ISFFT_PLAN_ESTIMATE) then
         call Estimate_Kernel(p)
         if (present(plan)) plan = p
         return
      end if
      select case (transform_kind)
      case (1)
         call ISFFT_Tune_FFT_R2C_Shape(transform_length, line_count)
      case (2)
         call ISFFT_Tune_FFT_Complex_Shape(transform_length, line_count)
      case (3, 4)
         p%transform_length = transform_length
         p%transform_kind = transform_kind
         p%line_count = line_count
         call Tune_Real_Kernel_And_Inner_Batch(p)
      end select
      if (transform_kind <= 2) then
         call Require_FFT_Precompute(transform_length, slot)
         if (transform_kind == 1) then
            do e = 1, fft_precompute_slots(slot)%r2c_strategy_count
               if (fft_precompute_slots(slot)%r2c_strategy_lines(e) /= line_count) cycle
               p%kernel_method = fft_precompute_slots(slot)%r2c_strategy_method(e)
               p%inner_batch_width = fft_precompute_slots(slot)%r2c_strategy_width(e)
            end do
         else
            do e = 1, fft_precompute_slots(slot)%complex_strategy_count
               if (fft_precompute_slots(slot)%complex_strategy_lines(e) /= line_count) cycle
               p%kernel_method = fft_precompute_slots(slot)%complex_strategy_method(e)
               p%inner_batch_width = fft_precompute_slots(slot)%complex_strategy_width(e)
            end do
         end if
      end if
      !! Partial MPI cache hits still execute collective tuning for missing ranks.
      !! Preserve this rank's previously returned numerical plan after participating.
      if (hit > 0) p = prepared_kernel_plans(hit)
      if (.not. allocated(prepared_kernel_plans)) allocate (prepared_kernel_plans(0))
      do e = 1, size(prepared_kernel_plans)
         if (prepared_kernel_plans(e)%transform_length == transform_length .and. &
            prepared_kernel_plans(e)%transform_kind == transform_kind .and. &
            prepared_kernel_plans(e)%line_count == line_count) exit
      end do
      if (e > size(prepared_kernel_plans)) then
         count = size(prepared_kernel_plans)
         allocate (grown(count + 1))
         grown(1:count) = prepared_kernel_plans
         call move_alloc(grown, prepared_kernel_plans)
      end if
      prepared_kernel_plans(e) = p
      if (present(plan)) plan = p
   end subroutine ! ISFFT_Prepare_Kernel_Impl

   subroutine ISFFT_Measure_FFT_Complex_Candidate( &
      f, transform_length, slot_index, &
      method, batch_width, repeat_count, global_time)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)
      integer, intent(in) :: transform_length, slot_index, method, batch_width, repeat_count
      real(kind=ISFFT_REAL_KIND), intent(out) :: global_time
      integer :: repeat, i, j, k, ierr
      real(kind=ISFFT_REAL_KIND) :: time_start, local_time

      do k = 1, size(f, 3)
         do j = 1, size(f, 2)
            do i = 1, transform_length
               f(i, j, k) = cmplx(real(mod(i + 3 * j + 5 * k, 97), ISFFT_REAL_KIND), &
                  real(mod(7 * i + j + 2 * k, 89), ISFFT_REAL_KIND), kind= &
                  ISFFT_REAL_KIND)
            end do
         end do
      end do
      call MPI_Barrier(ISFFT_Comm, ierr)
      time_start = MPI_Wtime()
      do repeat = 1, repeat_count
         call ISFFT_Execute_FFT_Complex_Strategy(transform_length, f, +1, slot_index, method, batch_width)
         call ISFFT_Execute_FFT_Complex_Strategy(transform_length, f, -1, slot_index, method, batch_width)
      end do
      local_time = MPI_Wtime() - time_start
      call MPI_Allreduce(local_time, global_time, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
   end subroutine ISFFT_Measure_FFT_Complex_Candidate

   subroutine ISFFT_Measure_FFT_R2C_Candidate( &
      full_re, spectrum, transform_length, &
      method, batch_width, repeat_count, global_time)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: full_re(:, :, :)
      complex(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: spectrum(:, :, :)
      integer, intent(in) :: transform_length, method, batch_width, repeat_count
      real(kind=ISFFT_REAL_KIND), intent(out) :: global_time
      integer :: repeat, i, j, k, ierr
      real(kind=ISFFT_REAL_KIND) :: time_start, local_time

      do k = 1, size(full_re, 3)
         do j = 1, size(full_re, 2)
            do i = 1, transform_length
               full_re(i, j, k) = real(mod(i + 3 * j + 5 * k, 97), ISFFT_REAL_KIND)
            end do
         end do
      end do
      spectrum = cmplx(0.0_ISFFT_REAL_KIND, 0.0_ISFFT_REAL_KIND, kind=ISFFT_REAL_KIND)
      call MPI_Barrier(ISFFT_Comm, ierr)
      time_start = MPI_Wtime()
      do repeat = 1, repeat_count
         if (method == 1) then
            call ISFFT_Execute_R2C_Scalar(transform_length, full_re, spectrum, +1)
            call ISFFT_Execute_R2C_Scalar(transform_length, full_re, spectrum, -1)
         else if (method == 5) then
            call ISFFT_Execute_R2C_Odd_Batch(transform_length, full_re, spectrum, +1, .false., batch_width)
            call ISFFT_Execute_R2C_Odd_Batch(transform_length, full_re, spectrum, -1, .false., batch_width)
         else if (method == 2 .or. method == 4) then
            call ISFFT_Execute_R2C_Batch(transform_length, full_re, spectrum, +1, batch_width)
            call ISFFT_Execute_R2C_Batch(transform_length, full_re, spectrum, -1, batch_width)
         else
            call ISFFT_Execute_R2C_Generic_Batch(transform_length, full_re, spectrum, +1, batch_width)
            call ISFFT_Execute_R2C_Generic_Batch(transform_length, full_re, spectrum, -1, batch_width)
         end if
      end do
      local_time = MPI_Wtime() - time_start
      call MPI_Allreduce(local_time, global_time, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
   end subroutine ISFFT_Measure_FFT_R2C_Candidate

   !! One numerical candidate executor for all tested grades. Full transform length
   !! is retained; only the number of independent lines is sampled.
   recursive subroutine Estimate_Kernel(plan)
      type(Kernel_Plan), intent(inout) :: plan
      type(Kernel_Plan) :: child_plan
      integer :: slot, child, limit, l

      plan%kernel_method = 1
      plan%inner_batch_width = 1
      if (plan%transform_kind == 2) then
         call Ensure_FFT_Precompute(plan%transform_length, slot)
         if (fft_precompute_slots(slot)%is_radix2) then
            plan%kernel_method = 4
         else if (fft_precompute_slots(slot)%use_rader) then
            !! Native scalar Rader owns its convolution; avoid unrelated mixed factors.
            plan%kernel_method = 1
         else
            plan%kernel_method = 5
         end if
         if (plan%kernel_method /= 1) then
            fft_precompute_slots(slot)%batch_ready = .true.
            call Reserve_FFT_Base(slot)
            plan%inner_batch_width = Estimate_Width(plan%transform_length, plan%line_count, 16, 3, 16)
         end if
         call ISFFT_Store_Complex_Strategy(slot, plan%line_count, plan%kernel_method, plan%inner_batch_width)
      else if (plan%transform_kind == 1) then
         call Ensure_R2C_Precompute(plan%transform_length, slot)
         if (mod(plan%transform_length, 2) == 0) then
            call Ensure_FFT_Precompute(plan%transform_length / 2, child)
            if (fft_precompute_slots(child)%is_radix2) then
               plan%kernel_method = 2
               plan%inner_batch_width = Estimate_Width(plan%transform_length, plan%line_count, 8, 4, 16)
            end if
         else
            plan%kernel_method = 5
            plan%inner_batch_width = Estimate_Width(plan%transform_length, max(1, (plan%line_count + 1) / 2), 16, 3, 64)
            child_plan%transform_length = plan%transform_length
            child_plan%transform_kind = 2
            child_plan%line_count = plan%inner_batch_width
            call Estimate_Kernel(child_plan)
            child_plan%line_count = max(1, mod((plan%line_count + 1) / 2, plan%inner_batch_width))
            call Estimate_Kernel(child_plan)
         end if
         call Reserve_R2C_Workspace(plan%transform_length, plan%line_count, plan%kernel_method)
         call ISFFT_Store_R2C_Strategy(slot, plan%line_count, plan%kernel_method, plan%inner_batch_width)
      else
         plan%kernel_method = 1
         plan%inner_batch_width = Estimate_Width(plan%transform_length, plan%line_count, 8, 4, 32)
         l = max(1, plan%transform_length - 1)
         call Prepare_Real_Workspace(l, plan%line_count)
         child_plan%transform_length = l
         child_plan%transform_kind = 2
         child_plan%line_count = (plan%inner_batch_width + 1) / 2
         if (mod(l, 2) == 0) then
            child_plan%transform_kind = 1
            child_plan%line_count = plan%inner_batch_width
         end if
         call Estimate_Kernel(child_plan)
         limit = max(1, mod(plan%line_count, plan%inner_batch_width))
         child_plan%line_count = limit
         if (mod(l, 2) /= 0) child_plan%line_count = (limit + 1) / 2
         call Estimate_Kernel(child_plan)
      end if
   end subroutine ! Estimate_Kernel

   recursive subroutine ISFFT_Tune_FFT_Complex_Shape(transform_length, line_count)
      integer, intent(in) :: transform_length, line_count

      if (transform_length <= 1) return
      if (Kernel_Shape_Cached(transform_length, 2, max(1, line_count))) return
      call Tune_FFT_Kernel_And_Inner_Batch(transform_length, 2, max(1, line_count))
   end subroutine ! ISFFT_Tune_FFT_Complex_Shape

   recursive subroutine ISFFT_Tune_FFT_R2C_Shape(transform_length, line_count)
      integer, intent(in) :: transform_length, line_count

      if (transform_length <= 1) return
      if (Kernel_Shape_Cached(transform_length, 1, max(1, line_count))) return
      call Tune_FFT_Kernel_And_Inner_Batch(transform_length, 1, max(1, line_count))
   end subroutine ! ISFFT_Tune_FFT_R2C_Shape

   subroutine ISFFT_Tune_FFT_Complex_Batch(f)
      complex(ISFFT_REAL_KIND), intent(inout), contiguous, target :: f(:, :, :)

      call ISFFT_Tune_FFT_Complex_Shape(size(f, 1), size(f, 2) * size(f, 3))
   end subroutine ! ISFFT_Tune_FFT_Complex_Batch

   subroutine ISFFT_Tune_FFT_R2C_Batch(r, z)
      real(ISFFT_REAL_KIND), intent(inout), contiguous, target :: r(:, :, :)
      complex(ISFFT_REAL_KIND), intent(inout), contiguous, target :: z(:, :, :)

      call ISFFT_Tune_FFT_R2C_Shape(size(r, 1), size(r, 2) * size(r, 3))
   end subroutine ! ISFFT_Tune_FFT_R2C_Batch

   !! Candidate method IDs must match the FFT_kernel execution dispatch.
   !! Keep ESTIMATE outside all timed candidate trials.
   recursive subroutine Tune_FFT_Kernel_And_Inner_Batch(transform_length, transform_kind, line_count)
      integer, intent(in) :: transform_length, transform_kind, line_count
      type(Kernel_Plan) :: base
      integer :: slot, child, ns, nc, c, m, w, j, rep, lim, anchor, global_lines, ierr, best, duplicate, &
         all_duplicate, baseline, &
         base_method, base_width
      integer :: methods(256), widths(256), reuse, k, local_hit(2), global_hit(2)
      type(Sample_Scores) :: saved_scores
      integer :: unique_local, unique_global
      integer, allocatable :: choices(:)
      real(ISFFT_REAL_KIND) :: scores(256), elapsed, total
      complex(ISFFT_REAL_KIND), allocatable :: z(:, :, :), sp(:, :, :)
      real(ISFFT_REAL_KIND), allocatable :: r(:, :, :)

      base%transform_length = transform_length
      base%transform_kind = transform_kind
      base%line_count = line_count
      if (planning_level == ISFFT_PLAN_ESTIMATE) then
         call Estimate_Kernel(base)
         return
      end if
      !! Estimate's starting width is included without storing its strategy as a tuned result.
      lim = 16
      if (transform_kind == 1 .and. mod(transform_length, 2) /= 0) lim = 64
      call MPI_Allreduce(line_count, global_lines, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      anchor = Estimate_Width(transform_length, global_lines, merge(16, 8, transform_kind == 2), &
         merge(3, 4, transform_kind == 2), lim)
      call Policy_Widths(planning_level, anchor, lim, merge(3, 1, lim == 64), choices)
      if (transform_kind == 2) then
         call Prepare_Complex_Candidates(transform_length, slot)
      else
         call Prepare_R2C_Candidates(transform_length, line_count, slot)
         if (mod(transform_length, 2) == 0) call Ensure_FFT_Precompute(transform_length / 2, child)
      end if
      nc = 0
      do m = 1, 12
         if (transform_kind == 2) then
            if ((m == 2 .or. m == 3) .and. .not. fft_precompute_slots(slot)%is_radix4) cycle
            if (m == 4 .and. .not. fft_precompute_slots(slot)%is_radix2) cycle
            if ((m == 5 .or. m == 11) .and. fft_precompute_slots(slot)%bluestein_m <= 0) cycle
            if (any(m == [6, 7, 9]) .and. fft_precompute_slots(slot)%is_radix2) cycle
            if ((m == 8 .or. m == 12) .and. .not. fft_precompute_slots(slot)%rader_ready) cycle
            if (m == 10 .and. .not. ISFFT_Is_Three_Times_Power_Of_Two(transform_length)) cycle
         else
            if (m > 5) cycle
            if (mod(transform_length, 2) /= 0) then
               if (m /= 1 .and. m /= 5) cycle
            else
               if (m == 5) cycle
               if (m == 2 .and. .not. fft_precompute_slots(child)%is_radix2) cycle
               if (m == 3 .and. .not. fft_precompute_slots(slot)%use_bluestein) cycle
               if (m == 4 .and. .not. ISFFT_Is_Three_Times_Power_Of_Two(transform_length / 2)) cycle
            end if
         end if
         do j = 1, size(choices)
            w = choices(j)
            if (m == 1 .and. w /= 1) cycle
            if (transform_kind == 2 .and. any(m == [2, 6, 11, 12]) .and. w /= 1) cycle
            nc = nc + 1
            methods(nc) = m
            widths(nc) = min(w, line_count)
            if (transform_kind == 1 .and. m == 5) widths(nc) = min(w, max(1, (line_count + 1) / 2))
         end do
      end do
      if (nc < 1) error stop 'no legal FFT tuning candidate'
      !! Selection is unnecessary when the legal candidate set has one element.
      if (planning_level == ISFFT_PLAN_MEASURE) then
         unique_local = 1
         do c = 2, nc
            if (methods(c) /= methods(1) .or. widths(c) /= widths(1)) unique_local = 0
         end do
         call MPI_Allreduce(unique_local, unique_global, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         if (unique_global == 1) then
            if (transform_kind == 2) then
               call ISFFT_Store_Complex_Strategy(slot, line_count, methods(1), widths(1))
            else
               !! Odd real candidates need their complex dependency before execution.
               if (mod(transform_length, 2) /= 0 .and. methods(1) == 5) then
                  call ISFFT_Tune_FFT_Complex_Shape(transform_length, widths(1))
                  call ISFFT_Tune_FFT_Complex_Shape(transform_length, max(1, mod((line_count + 1) / 2, widths(1))))
               end if
               call ISFFT_Store_R2C_Strategy(slot, line_count, methods(1), widths(1))
            end if
            return
         end if
      end if
      !! All ranks traverse the same candidate IDs, even when their local work differs.
      if (transform_kind == 1 .and. mod(transform_length, 2) /= 0) then
         do c = 1, nc
            if (methods(c) /= 5) cycle
            call ISFFT_Tune_FFT_Complex_Shape(transform_length, widths(c))
            call ISFFT_Tune_FFT_Complex_Shape(transform_length, max(1, mod((line_count + 1) / 2, widths(c))))
         end do
      end if
      ns = min(line_count, max(Policy_Sample(planning_level, line_count), 2 * maxval(widths(:nc))))
      reuse = 0
      if (planning_level == ISFFT_PLAN_MEASURE .and. transform_kind == 2) then
         if (allocated(sample_cache)) then
            do k = 1, size(sample_cache)
               associate (q => sample_cache(k))
                  if (q%transform_length /= transform_length .or. q%sample_line_count /= ns .or. &
                     q%candidate_count /= nc .or. q%planning_level /= planning_level) cycle
                  if (q%warmup_count /= Policy_Warmup_Count(planning_level, KERNEL_TRIAL) .or. &
                     q%repeat_count /= Policy_Repeat_Count(planning_level, KERNEL_TRIAL)) cycle
                  if (any(q%methods /= methods(:nc)) .or. any(q%widths /= widths(:nc))) cycle
                  reuse = k
                  exit
               end associate
            end do
         end if
         !! Require the same collective experiment ID, not unrelated local matches.
         local_hit = [reuse, -reuse]
         call MPI_Allreduce(local_hit, global_hit, 2, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         if (global_hit(1) <= 0 .or. global_hit(1) /= -global_hit(2)) reuse = 0
      end if
      scores = huge(1.0_ISFFT_REAL_KIND)
      if (reuse > 0) then
         scores(:nc) = sample_cache(reuse)%scores
      else
         if (transform_kind == 2) then
            allocate (z(transform_length, ns, 1))
         else
            allocate (r(transform_length, ns, 1), sp(transform_length / 2 + 1, ns, 1))
         end if
         do c = 1, nc
            duplicate = 0
            do j = 1, c - 1
               if (methods(c) == methods(j) .and. widths(c) == widths(j)) duplicate = 1
            end do
            call MPI_Allreduce(duplicate, all_duplicate, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
            if (all_duplicate == 1) cycle
            total = 0
            do rep = 1, Policy_Warmup_Count(planning_level, KERNEL_TRIAL) + &
               Policy_Repeat_Count(planning_level, KERNEL_TRIAL)
               if (transform_kind == 2) then
                  call ISFFT_Measure_FFT_Complex_Candidate( &
                     z, transform_length, slot, &
                     methods(c), widths(c), 1, elapsed)
               else
                  call ISFFT_Measure_FFT_R2C_Candidate( &
                     r, sp, transform_length, &
                     methods(c), widths(c), 1, elapsed)
               end if
               if (rep > Policy_Warmup_Count(planning_level, KERNEL_TRIAL)) total = total + elapsed
            end do
            scores(c) = total / Policy_Repeat_Count(planning_level, KERNEL_TRIAL)
         end do
         if (planning_level == ISFFT_PLAN_MEASURE .and. transform_kind == 2) then
            saved_scores%transform_length = transform_length
            saved_scores%sample_line_count = ns
            saved_scores%candidate_count = nc
            saved_scores%planning_level = planning_level
            saved_scores%warmup_count = Policy_Warmup_Count(planning_level, KERNEL_TRIAL)
            saved_scores%repeat_count = Policy_Repeat_Count(planning_level, KERNEL_TRIAL)
            saved_scores%methods = methods(:nc)
            saved_scores%widths = widths(:nc)
            saved_scores%scores = scores(:nc)
            if (.not. allocated(sample_cache)) allocate (sample_cache(0))
            sample_cache = [sample_cache, saved_scores]
         end if
      end if
      !! Reapply selection to the current execution workload, including its baseline.
      base_method = 1
      if (transform_kind == 2) then
         if (fft_precompute_slots(slot)%is_radix2) base_method = 4
         if (fft_precompute_slots(slot)%use_bluestein) base_method = 5
      else if (mod(transform_length, 2) /= 0) then
         base_method = 5
      else
         if (fft_precompute_slots(child)%is_radix2) base_method = 2
      end if
      base_width = Estimate_Width(transform_length, line_count, merge(16, 8, transform_kind == 2), &
         merge(3, 4, transform_kind == 2), lim)
      if (transform_kind == 1 .and. mod(transform_length, 2) /= 0) base_width = &
         Estimate_Width(transform_length, max(1, (line_count + 1) / 2), 16, 3, 64)
      baseline = minloc(scores(:nc), dim=1)
      do c = 1, nc
         if (methods(c) == base_method .and. widths(c) == base_width .and. scores(c) < huge(elapsed)) baseline = c
      end do
      best = Select_Best_Candidate(scores(:nc), baseline, real(MPI_Wtick(), ISFFT_REAL_KIND))
      if (transform_kind == 2) then
         call ISFFT_Store_Complex_Strategy(slot, line_count, methods(best), widths(best))
      else
         call ISFFT_Store_R2C_Strategy(slot, line_count, methods(best), widths(best))
      end if
   end subroutine ! Tune_FFT_Kernel_And_Inner_Batch

   subroutine Tune_Real_Kernel_And_Inner_Batch(plan)
      type(Kernel_Plan), intent(inout) :: plan
      integer, allocatable :: choices(:)
      integer :: global_lines, anchor, ierr, c, w, method, rep, ns, l, tail, best_width, best_method, i, j, &
         duplicate, all_duplicate
      integer, allocatable :: requests(:), skip_local(:), skip_global(:)
      integer :: d, q, request_lines, request_kind
      real(ISFFT_REAL_KIND), allocatable :: sample(:, :, :), v(:, :, :)
      real(ISFFT_REAL_KIND) :: t, dt, elapsed, score, best

      call MPI_Allreduce(plan%line_count, global_lines, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      anchor = Estimate_Width(plan%transform_length, global_lines, 8, 4, 32)
      call Policy_Widths(planning_level, anchor, 64, 2, choices)
      choices = [0, choices]
      l = max(1, plan%transform_length - 1)
      if (planning_level == ISFFT_PLAN_MEASURE) then
         allocate (requests(2 * size(choices)), skip_local(2 * size(choices)), skip_global(2 * size(choices)))
         do c = 1, size(choices)
            w = max(1, min(choices(c), plan%line_count))
            tail = max(1, mod(plan%line_count, w))
            do d = 1, 2
               request_lines = w
               if (d == 2) request_lines = tail
               if (mod(l, 2) /= 0) request_lines = (request_lines + 1) / 2
               q = 2 * (c - 1) + d
               requests(q) = request_lines
               skip_local(q) = 0
               if (any(requests(:q - 1) == request_lines)) skip_local(q) = 1
            end do
         end do
         !! All ranks agree once. A request is skipped only if every rank has
         !! already requested its own identical dependency in this traversal.
         call MPI_Allreduce(skip_local, skip_global, size(requests), MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      end if
      do c = 1, size(choices)
         w = max(1, min(choices(c), plan%line_count))
         tail = max(1, mod(plan%line_count, w))
         do d = 1, 2
            request_lines = w
            if (d == 2) request_lines = tail
            request_kind = 1
            if (mod(l, 2) /= 0) then
               request_kind = 2
               request_lines = (request_lines + 1) / 2
            end if
            if (planning_level == ISFFT_PLAN_MEASURE) then
               if (skip_global(2 * (c - 1) + d) == 1) then
                  cycle
               end if
            end if
            if (request_kind == 1) then
               call ISFFT_Tune_FFT_R2C_Shape(l, request_lines)
            else
               call ISFFT_Tune_FFT_Complex_Shape(l, request_lines)
            end if
         end do
      end do
      ns = min(plan%line_count, max(Policy_Sample(planning_level, plan%line_count), 64))
      allocate (sample(plan%transform_length, ns, 1), v(plan%transform_length, ns, 1))
      do j = 1, ns
         do i = 1, plan%transform_length
            sample(i, j, 1) = sin(real(i + 7 * j, ISFFT_REAL_KIND))
         end do
      end do
      if (plan%transform_kind == 3) then
         sample(1, :, :) = 0
         sample(plan%transform_length, :, :) = 0
      end if
      best = huge(best)
      best_method = 0
      best_width = 1
      do c = 1, size(choices)
         w = max(1, min(choices(c), plan%line_count))
         method = 1
         if (choices(c) == 0) method = 0
         duplicate = 0
         if (c > 2) then
            if (any(min(choices(2:c - 1), plan%line_count) == w)) duplicate = 1
         end if
         call MPI_Allreduce(duplicate, all_duplicate, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         if (all_duplicate == 1) cycle
         score = 0
         do rep = 1, Policy_Warmup_Count(planning_level, KERNEL_TRIAL) + &
            Policy_Repeat_Count(planning_level, KERNEL_TRIAL)
            v = sample
            call MPI_Barrier(ISFFT_Comm, ierr)
            t = MPI_Wtime()
            if (plan%transform_kind == 3) then
               call ISFFT_Execute_DST_Plan(plan%transform_length, v, 1, method, w)
               call ISFFT_Execute_DST_Plan(plan%transform_length, v, -1, method, w)
            else
               call ISFFT_Execute_DCT_Plan(plan%transform_length, v, 1, method, w)
               call ISFFT_Execute_DCT_Plan(plan%transform_length, v, -1, method, w)
            end if
            dt = MPI_Wtime() - t
            call MPI_Allreduce(dt, elapsed, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
            if (rep > Policy_Warmup_Count(planning_level, KERNEL_TRIAL)) score = score + elapsed
         end do
         if (score < best) then
            best = score
            best_width = w
            best_method = method
         end if
      end do
      plan%inner_batch_width = best_width
      plan%kernel_method = best_method
   end subroutine ! Tune_Real_Kernel_And_Inner_Batch

   subroutine ISFFT_Store_Complex_Strategy(slot_index, line_count, method, width)
      integer, intent(in) :: slot_index, line_count, method, width
      integer, allocatable :: lines(:), methods(:), widths(:)
      integer :: old_count, new_count, entry

      old_count = fft_precompute_slots(slot_index)%complex_strategy_count
      !! Keep a previously returned choice stable on partial MPI cache hits.
      do entry = 1, old_count
         if (fft_precompute_slots(slot_index)%complex_strategy_lines(entry) == line_count) return
      end do
      new_count = old_count + 1
      allocate (lines(new_count), methods(new_count), widths(new_count))
      if (old_count > 0) then
         lines(1:old_count) = fft_precompute_slots(slot_index)%complex_strategy_lines(1:old_count)
         methods(1:old_count) = fft_precompute_slots(slot_index)%complex_strategy_method(1:old_count)
         widths(1:old_count) = fft_precompute_slots(slot_index)%complex_strategy_width(1:old_count)
      end if
      lines(new_count) = line_count
      methods(new_count) = method
      widths(new_count) = width
      call move_alloc(lines, fft_precompute_slots(slot_index)%complex_strategy_lines)
      call move_alloc(methods, fft_precompute_slots(slot_index)%complex_strategy_method)
      call move_alloc(widths, fft_precompute_slots(slot_index)%complex_strategy_width)
      fft_precompute_slots(slot_index)%complex_strategy_count = new_count
   end subroutine ISFFT_Store_Complex_Strategy

   subroutine ISFFT_Store_R2C_Strategy(slot_index, line_count, method, width)
      integer, intent(in) :: slot_index, line_count, method, width
      integer, allocatable :: lines(:), methods(:), widths(:)
      integer :: old_count, new_count, entry

      old_count = fft_precompute_slots(slot_index)%r2c_strategy_count
      !! Keep a previously returned choice stable on partial MPI cache hits.
      do entry = 1, old_count
         if (fft_precompute_slots(slot_index)%r2c_strategy_lines(entry) == line_count) return
      end do
      new_count = old_count + 1
      allocate (lines(new_count), methods(new_count), widths(new_count))
      if (old_count > 0) then
         lines(1:old_count) = fft_precompute_slots(slot_index)%r2c_strategy_lines(1:old_count)
         methods(1:old_count) = fft_precompute_slots(slot_index)%r2c_strategy_method(1:old_count)
         widths(1:old_count) = fft_precompute_slots(slot_index)%r2c_strategy_width(1:old_count)
      end if
      lines(new_count) = line_count
      methods(new_count) = method
      widths(new_count) = width
      call move_alloc(lines, fft_precompute_slots(slot_index)%r2c_strategy_lines)
      call move_alloc(methods, fft_precompute_slots(slot_index)%r2c_strategy_method)
      call move_alloc(widths, fft_precompute_slots(slot_index)%r2c_strategy_width)
      fft_precompute_slots(slot_index)%r2c_strategy_count = new_count
   end subroutine ISFFT_Store_R2C_Strategy

   subroutine ISFFT_Clone_Complex_Strategy(transform_length, source_lines, target_lines)
      integer, intent(in) :: transform_length, source_lines, target_lines
      integer :: slot_n, entry, source_entry

      call Require_FFT_Precompute(transform_length, slot_n)
      source_entry = 0
      do entry = 1, fft_precompute_slots(slot_n)%complex_strategy_count
         if (fft_precompute_slots(slot_n)%complex_strategy_lines(entry) == target_lines) return
         if (fft_precompute_slots(slot_n)%complex_strategy_lines(entry) == source_lines) source_entry = entry
      end do
      if (source_entry == 0) error stop 'ISFFT_Clone_Complex_Strategy: source plan missing'
      call ISFFT_Store_Complex_Strategy(slot_n, target_lines, &
         fft_precompute_slots(slot_n)%complex_strategy_method(source_entry), &
         fft_precompute_slots(slot_n)%complex_strategy_width(source_entry))
   end subroutine ISFFT_Clone_Complex_Strategy

   subroutine ISFFT_Clone_R2C_Strategy(transform_length, source_lines, target_lines)
      integer, intent(in) :: transform_length, source_lines, target_lines
      integer :: slot_n, entry, source_entry

      call Require_FFT_Precompute(transform_length, slot_n)
      source_entry = 0
      do entry = 1, fft_precompute_slots(slot_n)%r2c_strategy_count
         if (fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) == target_lines) return
         if (fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) == source_lines) source_entry = entry
      end do
      if (source_entry == 0) error stop 'ISFFT_Clone_R2C_Strategy: source plan missing'
      call ISFFT_Store_R2C_Strategy(slot_n, target_lines, &
         fft_precompute_slots(slot_n)%r2c_strategy_method(source_entry), &
         fft_precompute_slots(slot_n)%r2c_strategy_width(source_entry))
   end subroutine ISFFT_Clone_R2C_Strategy

   subroutine ISFFT_Forget_Complex_Strategy(transform_length, line_count)
      integer, intent(in) :: transform_length, line_count
      integer :: slot_n, entry, last

      call Require_FFT_Precompute(transform_length, slot_n)
      last = fft_precompute_slots(slot_n)%complex_strategy_count
      do entry = 1, last
         if (fft_precompute_slots(slot_n)%complex_strategy_lines(entry) /= line_count) cycle
         fft_precompute_slots(slot_n)%complex_strategy_lines(entry) = &
            fft_precompute_slots(slot_n)%complex_strategy_lines(last)
         fft_precompute_slots(slot_n)%complex_strategy_method(entry) = &
            fft_precompute_slots(slot_n)%complex_strategy_method(last)
         fft_precompute_slots(slot_n)%complex_strategy_width(entry) = &
            fft_precompute_slots(slot_n)%complex_strategy_width(last)
         fft_precompute_slots(slot_n)%complex_strategy_count = last - 1
         return
      end do
   end subroutine ISFFT_Forget_Complex_Strategy

   subroutine ISFFT_Forget_R2C_Strategy(transform_length, line_count)
      integer, intent(in) :: transform_length, line_count
      integer :: slot_n, entry, last

      call Require_FFT_Precompute(transform_length, slot_n)
      last = fft_precompute_slots(slot_n)%r2c_strategy_count
      do entry = 1, last
         if (fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) /= line_count) cycle
         fft_precompute_slots(slot_n)%r2c_strategy_lines(entry) = fft_precompute_slots(slot_n)%r2c_strategy_lines(last)
         fft_precompute_slots(slot_n)%r2c_strategy_method(entry) = &
            fft_precompute_slots(slot_n)%r2c_strategy_method(last)
         fft_precompute_slots(slot_n)%r2c_strategy_width(entry) = fft_precompute_slots(slot_n)%r2c_strategy_width(last)
         fft_precompute_slots(slot_n)%r2c_strategy_count = last - 1
         return
      end do
   end subroutine ISFFT_Forget_R2C_Strategy

   logical function Kernel_Shape_Cached(transform_length, transform_kind, line_count) result(cached)
      integer, intent(in) :: transform_length, transform_kind, line_count
      integer :: slot, entry, hit, all_hit, ierr

      if (transform_kind == 1) then
         call Ensure_R2C_Precompute(transform_length, slot)
      else
         call Ensure_FFT_Precompute(transform_length, slot)
      end if
      hit = 0
      if (transform_kind == 1) then
         do entry = 1, fft_precompute_slots(slot)%r2c_strategy_count
            if (fft_precompute_slots(slot)%r2c_strategy_lines(entry) == line_count) hit = 1
         end do
      else
         do entry = 1, fft_precompute_slots(slot)%complex_strategy_count
            if (fft_precompute_slots(slot)%complex_strategy_lines(entry) == line_count) hit = 1
         end do
      end if
      call MPI_Allreduce(hit, all_hit, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      cached = all_hit == 1
   end function ! Kernel_Shape_Cached

   !------------------------------------------------------------------------------
   !! Access, outer-batch and distributed-block selection
   !------------------------------------------------------------------------------
   subroutine Plan_Transport(plan, choice)
      type(Transfer_Plan), intent(inout) :: plan
      integer, intent(in) :: choice
      integer :: selected

      selected = transport_policy
      if (choice == 1 .and. selected < 0) selected = 0
      if (choice == 2 .and. selected < 0) selected = 1
      if (selected < 0) then
         call Plan_Communication_Strategy(plan)
      else
         !! Public policy uses the same IDs as Select_Transfer_Method.
         call Select_Transfer_Method(plan, selected)
      end if
   end subroutine ! Plan_Transport

   subroutine Tune_Block_Width(plan)
      type(Axis_Plan), intent(inout) :: plan
      type(Block_Choice), allocatable :: grown(:)
      integer :: i, hit, allhit, ierr, count, method, width, lines

      if (planning_level == ISFFT_PLAN_ESTIMATE) then
         call Tune_Block_Width_Impl(plan)
         return
      end if
      hit = 0
      method = plan%kernel%kernel_method
      width = plan%kernel%inner_batch_width
      lines = plan%kernel%line_count
      if (allocated(block_cache)) then
         do i = 1, size(block_cache)
            associate (q => block_cache(i)%plan)
               if (plan%physical_axis /= q%physical_axis .or. plan%transform_length /= q%transform_length .or. &
                  plan%transform_kind /= q%transform_kind) cycle
               if (method /= block_cache(i)%method .or. width /= block_cache(i)%width) cycle
               if (lines /= block_cache(i)%line_count) cycle
               if (.not. Same_Layout(plan%before%source, q%before%source) .or. &
                  .not. Same_Layout(plan%after%destination, q%after%destination)) cycle
               if (Same_Layout(plan%input, q%input) .and. Same_Layout(plan%output, q%output)) hit = i
            end associate
         end do
      end if
      call MPI_Allreduce(hit, allhit, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      if (allhit > 0) then
         plan%kernel = block_cache(hit)%plan%kernel
         plan%tail_kernel = block_cache(hit)%plan%tail_kernel
         plan%communication_batch_width = block_cache(hit)%plan%communication_batch_width
         return
      end if
      call Tune_Block_Width_Impl(plan)
      if (.not. allocated(block_cache)) allocate (block_cache(0))
      count = size(block_cache)
      allocate (grown(count + 1))
      grown(:count) = block_cache
      grown(count + 1)%plan = plan
      grown(count + 1)%method = method
      grown(count + 1)%width = width
      grown(count + 1)%line_count = lines
      call move_alloc(grown, block_cache)
   end subroutine ! Tune_Block_Width

   !! Estimate selects a block width from capacity and working-set rules.
   !! Measure prepares the kernel and tail plan associated with each candidate,
   !! then times the corresponding distributed block execution.
   !! Subsequent solves reuse the selected execution plans.
   subroutine Tune_Block_Width_Impl(plan)
      type(Axis_Plan), intent(inout) :: plan
      type(Kernel_Plan), allocatable :: kernels(:), tails(:)
      type(Spectral_Field) :: f
      real(ISFFT_REAL_KIND), allocatable :: original(:, :, :), scores(:)
      integer, allocatable :: widths(:)
      integer :: rank, ierr, n(3), off(3), lines, limit, global_limit, parts, bytes, anchor, c, best, i, j, k
      integer :: budget, warm, repeats

      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call Layout_Box(plan%before%source, rank, n, off)
      lines = product(pack(n, [1, 2, 3] /= plan%physical_axis))
      limit = Policy_Block_Capacity(plan%transform_length, lines)
      call MPI_Comm_size(plan%before%gather_comm, parts, ierr)
      bytes = merge(16, 8, plan%input%complex_data .or. plan%transform_kind == 1)
      if (planning_level == ISFFT_PLAN_ESTIMATE) then
         plan%communication_batch_width = Policy_Block_Start(plan%transform_length, max(1, lines), bytes, parts, limit)
         call Prepare_Block_Kernels(plan, lines)
         return
      end if
      call MPI_Allreduce(limit, global_limit, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      anchor = Policy_Block_Start(plan%transform_length, max(1, lines), bytes, parts, limit)
      call MPI_Allreduce(MPI_IN_PLACE, anchor, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      call Policy_Widths(planning_level, anchor, global_limit, 4, widths)
      allocate (kernels(size(widths)), tails(size(widths)), scores(size(widths)))
      !! Prepare dependencies before timed candidate execution. Existing exact-shape
      !! and sample-score caches share compatible requests within this initialization.
      do c = 1, size(widths)
         plan%communication_batch_width = min(widths(c), limit)
         call Prepare_Block_Kernels(plan, lines)
         kernels(c) = plan%kernel
         tails(c) = plan%tail_kernel
         call Reserve_Block(plan)
      end do
      allocate (original(n(1), n(2), n(3)))
      do k = 1, n(3)
         do j = 1, n(2)
            do i = 1, n(1)
               original(i, j, k) = sin(real(i + 7 * j + 13 * k + rank, ISFFT_REAL_KIND))
            end do
         end do
      end do
      if (plan%transform_kind == 3) then
         select case (plan%physical_axis)
         case (1)
            if (off(1) == 0 .and. n(1) > 0) original(1, :, :) = 0
            if (off(1) + n(1) == plan%transform_length .and. n(1) > 0) original(n(1), :, :) = 0
         case (2)
            if (off(2) == 0 .and. n(2) > 0) original(:, 1, :) = 0
            if (off(2) + n(2) == plan%transform_length .and. n(2) > 0) original(:, n(2), :) = 0
         case (3)
            if (off(3) == 0 .and. n(3) > 0) original(:, :, 1) = 0
            if (off(3) + n(3) == plan%transform_length .and. n(3) > 0) original(:, :, n(3)) = 0
         end select
      end if
      if (plan%before%source%complex_data) then
         allocate (f%z(size(original, 1), size(original, 2), size(original, 3)), &
            source=cmplx(original, original / 2, ISFFT_REAL_KIND))
      else
         allocate (f%r(size(original, 1), size(original, 2), size(original, 3)), source=original)
      end if
      !! Both directions reuse allocated fields even on the first measured trial.
      call Reserve_Field(f, plan%before%source)
      call Reserve_Field(f, plan%after%destination)
      budget = max(256, 2 * global_limit)
      warm = Policy_Warmup_Count(planning_level, SAMPLED_BATCH_TRIAL)
      repeats = Policy_Repeat_Count(planning_level, SAMPLED_BATCH_TRIAL)
      do c = 1, size(widths)
         call select_candidate(c)
         call evaluate(budget, SAMPLED_BATCH_TRIAL, scores(c))
      end do
      best = minloc(scores, dim=1)
      call select_candidate(best)
      call Clear_Field(f)
   contains

      subroutine select_candidate(index)
         integer, intent(in) :: index

         plan%communication_batch_width = min(widths(index), limit)
         plan%kernel = kernels(index)
         plan%tail_kernel = tails(index)
      end subroutine ! select_candidate

      subroutine evaluate(sample_budget, stage, mean_score)
         integer, intent(in) :: sample_budget, stage
         real(ISFFT_REAL_KIND), intent(out) :: mean_score
         integer :: rep, sampled
         real(ISFFT_REAL_KIND) :: t, elapsed, raw, estimated

         sampled = Block_Sample_Lines(lines, plan%communication_batch_width, sample_budget)
         mean_score = 0
         do rep = 1, warm + repeats
            if (plan%before%source%complex_data) then
               f%z = cmplx(original, original / 2, ISFFT_REAL_KIND)
            else
               f%r = original
            end if
            call MPI_Barrier(ISFFT_Comm, ierr)
            t = MPI_Wtime()
            call Execute_Block_Axis(plan, f, 1, sample_budget)
            call Execute_Block_Axis(plan, f, -1, sample_budget)
            elapsed = MPI_Wtime() - t
            !! Different batch widths can round the sampled line count differently.
            !! Compare projected local work, then take the slowest process, not raw
            !! times for unequal amounts of work. Full review has a factor of one.
            estimated = elapsed * real(max(1, lines), ISFFT_REAL_KIND) / max(1, sampled)
            call MPI_Allreduce(elapsed, raw, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
            call MPI_Allreduce(estimated, t, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
            if (rep > warm) mean_score = mean_score + t
         end do
         mean_score = mean_score / repeats
      end subroutine ! evaluate
   end subroutine ! Tune_Block_Width_Impl

   subroutine Prepare_Block_Kernels(plan, line_count)
      type(Axis_Plan), intent(inout) :: plan
      integer, intent(in) :: line_count
      integer :: width

      width = plan%communication_batch_width
      call Build_Kernel_Plan(plan%transform_length, plan%transform_kind, min(width, max(1, line_count)), plan%kernel, &
         plan%input%complex_data)
      call Tune_Access_And_Batching(plan%kernel, [plan%transform_length, plan%kernel%line_count, 1], 1)
      call Prepare_Block_Tail(plan, line_count, width)
   end subroutine ! Prepare_Block_Kernels

   subroutine Prepare_Block_Tail(plan, line_count, width)
      type(Axis_Plan), intent(inout) :: plan
      integer, intent(in) :: line_count, width
      integer :: tail, any_tail, ierr

      tail = mod(line_count, width)
      call MPI_Allreduce(tail, any_tail, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      plan%tail_kernel = Kernel_Plan()
      if (any_tail == 0) return
      !! Empty-tail ranks must still join collective kernel tuning when another
      !! rank has a real tail; only the all-rank absent case can skip the work.
      call Build_Kernel_Plan(plan%transform_length, plan%transform_kind, max(1, tail), plan%tail_kernel, &
         plan%input%complex_data)
      call Tune_Access_And_Batching(plan%tail_kernel, [plan%transform_length, plan%tail_kernel%line_count, 1], 1)
   end subroutine ! Prepare_Block_Tail

   subroutine Build_Kernel_Plan(transform_length, transform_kind, line_count, plan, complex_input)
      integer, intent(in) :: transform_length, transform_kind, line_count
      type(Kernel_Plan), intent(out) :: plan
      logical, intent(in), optional :: complex_input

      !! Preparation returns the numerical choice; Planner owns the outer-axis choice.
      call ISFFT_Prepare_Kernel(transform_length, transform_kind, max(1, line_count), plan=plan, axis_storage=.false.)
      if (present(complex_input)) plan%complex_input = complex_input

   end subroutine ! Build_Kernel_Plan

   !! Compare scalar vs genuine DST/DCT batches, independently of outer grouping.
   !! Outer grouping compares complete gather/call/scatter, with the chosen kernel fixed.
   subroutine Tune_Access_And_Batching(plan, local_shape, memory_axis)
      type(Kernel_Plan), intent(inout) :: plan
      integer, intent(in) :: local_shape(3), memory_axis

      call Tune_Access_And_Batching_Impl(plan, local_shape, memory_axis)
   end subroutine ! Tune_Access_And_Batching

   subroutine Tune_Access_And_Batching_Impl(plan, local_shape, memory_axis)
      type(Kernel_Plan), intent(inout) :: plan
      integer, intent(in) :: local_shape(3), memory_axis
      integer :: lines, global_lines, e, c, ierr, hit, all_hit, anchor, qshape(3), qlines, full_sample, all_full
      integer :: shortlist(4), which, nshort, best_c, cache_count, duplicate, all_duplicate, j
      integer, allocatable :: widths(:), selected(:)
      integer :: direct_anchor, stage, idx, pos, nsample
      logical :: anchor_direct_access, is_structural, is_full_batch
      type(Access_Candidate) :: structural_anchor
      type(Access_Candidate), allocatable :: candidates(:)
      real(ISFFT_REAL_KIND), allocatable :: scores(:), r(:, :, :), sample(:, :, :)
      complex(ISFFT_REAL_KIND), allocatable :: z(:, :, :), zs(:, :, :), spectrum(:, :, :)
      type(Kernel_Plan), allocatable :: grown(:)

      plan%local_shape = local_shape
      plan%memory_axis = memory_axis
      lines = int(product(int(local_shape, ISFFT_INT64_KIND)) / plan%transform_length)
      call Policy_Axis_Start(planning_level, plan%transform_length, plan%transform_kind, local_shape, memory_axis, &
         plan%complex_input, anchor, anchor_direct_access, direct_anchor)
      plan%outer_batch_width = anchor
      plan%direct_access = anchor_direct_access
      plan%direct_batch_width = direct_anchor
      if (planning_level == ISFFT_PLAN_ESTIMATE) return
      hit = 0
      if (allocated(axis_plan_cache)) then
         do e = 1, size(axis_plan_cache)
            associate (q => axis_plan_cache(e))
               if (q%transform_length == plan%transform_length .and. q%transform_kind == plan%transform_kind .and. &
                  q%kernel_method == plan%kernel_method .and. &
                  q%inner_batch_width == plan%inner_batch_width .and. all(q%local_shape == local_shape) .and. &
                  q%memory_axis == memory_axis .and. &
                  (q%complex_input .eqv. plan%complex_input)) hit = e
            end associate
         end do
      end if
      e = 0
      if (hit > 0) e = 1
      call MPI_Allreduce(e, all_hit, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      if (all_hit == 1) then
         plan%outer_batch_width = axis_plan_cache(hit)%outer_batch_width
         plan%direct_access = axis_plan_cache(hit)%direct_access
         plan%direct_batch_width = axis_plan_cache(hit)%direct_batch_width
         return
      end if
      call MPI_Allreduce(max(1, lines), global_lines, 1, MPI_INTEGER, MPI_MAX, ISFFT_Comm, ierr)
      anchor = Estimate_Width(plan%transform_length, global_lines, &
         merge(16, 8, plan%transform_kind == 2 .or. plan%complex_input), 3, global_lines, memory_axis)
      if (memory_axis == 1) anchor = global_lines
      if (planning_level == 1) then
         call Policy_Axis_Start(planning_level, plan%transform_length, plan%transform_kind, local_shape, memory_axis, &
            plan%complex_input, anchor, anchor_direct_access, direct_anchor)
         call MPI_Allreduce(MPI_IN_PLACE, anchor, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         call MPI_Allreduce(MPI_IN_PLACE, direct_anchor, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      end if
      call Policy_Axis_Widths(planning_level, anchor, global_lines, widths, plan%transform_length, &
         merge(16, 8, plan%transform_kind == 2 .or. plan%complex_input), &
         memory_axis)
      candidates = [(Access_Candidate(.false., widths(c)), c = 1, size(widths))]
      if (plan%transform_kind == 2 .and. memory_axis /= 1) &
         candidates = [candidates, Access_Candidate(.true., 1), Access_Candidate(.true., 4), &
         Access_Candidate(.true., 16)]
      structural_anchor = Access_Candidate(.false., anchor)
      if (planning_level == 1 .and. anchor_direct_access) then
         structural_anchor = Access_Candidate(.true., direct_anchor)
         if (.not. any(candidates%direct_access .and. candidates%batch_width == direct_anchor)) &
            candidates = [candidates, structural_anchor]
      end if
      allocate (scores(size(candidates)))
      scores = huge(1.0_ISFFT_REAL_KIND)
      stage = SAMPLED_BATCH_TRIAL
      qshape = local_shape
      if (memory_axis == 1 .and. lines > 0) qshape = &
         [plan%transform_length, Policy_Sample(planning_level, lines, .true.), 1]
      if (memory_axis == 2 .and. local_shape(1) > 0) qshape(3) = &
         min(local_shape(3), max(1, Policy_Sample(planning_level, lines, .true.) / local_shape(1)))
      qlines = int(product(int(qshape, ISFFT_INT64_KIND)) / plan%transform_length)
      if (planning_level == 1 .and. memory_axis == 3 .and. (plan%transform_kind == 2 .or. plan%complex_input) .and. &
         lines > 0) then
         nsample = min(lines, Policy_Sample(planning_level, lines, .true.))
         allocate (selected(nsample))
         !! Four equally spaced contiguous strips retain real strides and include the tail.
         do idx = 1, nsample
            pos = (idx - 1) / max(1, (nsample + 3) / 4)
            selected(idx) = 1 + min(pos, 3) * max(0, lines - nsample) / 3 + idx - 1
         end do
         qlines = nsample
      end if
      full_sample = merge(1, 0, qlines == lines)
      call MPI_Allreduce(full_sample, all_full, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
      call prepare_sample(qshape)
      do c = 1, size(candidates)
         duplicate = 0
         do j = 1, c - 1
            if (.not. candidates(c)%direct_access .and. .not. candidates(j)%direct_access) then
               if (min(candidates(c)%batch_width, qlines) == min(candidates(j)%batch_width, qlines)) duplicate = 1
            end if
         end do
         call MPI_Allreduce(duplicate, all_duplicate, 1, MPI_INTEGER, MPI_MIN, ISFFT_Comm, ierr)
         if (all_duplicate == 1) cycle
         call evaluate(candidates(c), scores(c))
      end do
      if (all_full == 0) then
         stage = FULL_BATCH_TRIAL
         nshort = 0
         do which = 1, 2
            best_c = minloc(scores, dim=1)
            if (scores(best_c) == huge(1.0_ISFFT_REAL_KIND)) exit
            nshort = nshort + 1
            shortlist(nshort) = best_c
            scores(best_c) = huge(1.0_ISFFT_REAL_KIND)
         end do
         !! Retain structural anchors even when a small sample aliases wide submissions.
         do c = 1, size(candidates)
            is_structural = (candidates(c)%direct_access .eqv. structural_anchor%direct_access) .and. &
               candidates(c)%batch_width == structural_anchor%batch_width
            is_full_batch = .not. candidates(c)%direct_access .and. candidates(c)%batch_width == global_lines
            if (.not. (is_structural .or. is_full_batch)) cycle
            if (any(shortlist(:nshort) == c)) cycle
            nshort = nshort + 1
            shortlist(nshort) = c
         end do
         scores = huge(1.0_ISFFT_REAL_KIND)
         if (allocated(selected)) deallocate (selected)
         qshape = local_shape
         qlines = lines
         call prepare_sample(local_shape)
         do which = 1, nshort
            c = shortlist(which)
            call evaluate(candidates(c), scores(c))
         end do
      end if
      best_c = minloc(scores, dim=1)
      call select_access(candidates(best_c))
      if (hit > 0) then
         plan%outer_batch_width = axis_plan_cache(hit)%outer_batch_width
         plan%direct_access = axis_plan_cache(hit)%direct_access
         plan%direct_batch_width = axis_plan_cache(hit)%direct_batch_width
      else
         if (.not. allocated(axis_plan_cache)) allocate (axis_plan_cache(0))
         cache_count = size(axis_plan_cache)
         allocate (grown(cache_count + 1))
         grown(:cache_count) = axis_plan_cache
         grown(cache_count + 1) = plan
         call move_alloc(grown, axis_plan_cache)
      end if
   contains

      subroutine prepare_sample(sz)
         integer, intent(in) :: sz(3)
         integer :: i, j, k

         if (allocated(sample)) deallocate (sample)
         if (allocated(r)) deallocate (r)
         if (allocated(z)) deallocate (z, zs)
         if (allocated(spectrum)) deallocate (spectrum)
         allocate (sample(sz(1), sz(2), sz(3)))
         do k = 1, sz(3)
            do j = 1, sz(2)
               do i = 1, sz(1)
                  sample(i, j, k) = sin(real(i + 3 * j + 7 * k, ISFFT_REAL_KIND))
               end do
            end do
         end do
         if (plan%transform_kind == 3 .and. size(sample) > 0) then
            select case (memory_axis)
            case (1)
               sample(1, :, :) = 0
               sample(plan%transform_length, :, :) = 0
            case (2)
               sample(:, 1, :) = 0
               sample(:, plan%transform_length, :) = 0
            case (3)
               sample(:, :, 1) = 0
               sample(:, :, plan%transform_length) = 0
            end select
         end if
         if (plan%transform_kind == 2 .or. (plan%transform_kind >= 3 .and. plan%complex_input)) then
            allocate (zs(size(sample, 1), size(sample, 2), size(sample, 3)), &
               source=cmplx(sample, 0.5_ISFFT_REAL_KIND * sample, ISFFT_REAL_KIND))
            allocate (z(size(zs, 1), size(zs, 2), size(zs, 3)), source=zs)
            deallocate (sample)
         else
            allocate (r(size(sample, 1), size(sample, 2), size(sample, 3)), source=sample)
         end if
         if (plan%transform_kind == 1) allocate (spectrum(plan%transform_length / 2 + 1, sz(2), sz(3)))
      end subroutine ! prepare_sample

      subroutine select_access(candidate)
         type(Access_Candidate), intent(in) :: candidate

         plan%direct_access = candidate%direct_access
         plan%outer_batch_width = 1
         plan%direct_batch_width = 1
         if (candidate%direct_access) then
            plan%direct_batch_width = candidate%batch_width
         else
            plan%outer_batch_width = max(1, min(max(1, lines), candidate%batch_width))
         end if
      end subroutine ! select_access

      subroutine evaluate(candidate, elapsed)
         type(Access_Candidate), intent(in) :: candidate
         real(ISFFT_REAL_KIND), intent(out) :: elapsed
         integer :: pass
         real(ISFFT_REAL_KIND) :: t, dt, gdt

         call select_access(candidate)
         call Reserve_Axis_Scratch(plan)
         elapsed = 0
         do pass = 1, Policy_Warmup_Count(planning_level, stage) + Policy_Repeat_Count(planning_level, stage)
            if (allocated(r)) r = sample
            if (allocated(z)) z = zs
            call MPI_Barrier(ISFFT_Comm, ierr)
            t = MPI_Wtime()
            select case (plan%transform_kind)
            case (1)
               call Execute_Kernel(plan, r, spectrum, 1)
               call Execute_Kernel(plan, r, spectrum, -1)
            case (2)
               call Execute_Kernel(plan, z, 1, memory_axis, selected)
               call Execute_Kernel(plan, z, -1, memory_axis, selected)
            case default
               if (plan%complex_input) then
                  call Execute_Kernel(plan, z, 1, memory_axis, selected)
                  call Execute_Kernel(plan, z, -1, memory_axis, selected)
               else
                  call Execute_Kernel(plan, r, 1)
                  call Execute_Kernel(plan, r, -1)
               end if
            end select
            dt = MPI_Wtime() - t
            call MPI_Allreduce(dt, gdt, 1, ISFFT_REAL_MPI_TYPE, MPI_MAX, ISFFT_Comm, ierr)
            if (pass > Policy_Warmup_Count(planning_level, stage)) elapsed = elapsed + gdt
         end do
         elapsed = elapsed / Policy_Repeat_Count(planning_level, stage)
      end subroutine ! evaluate
   end subroutine ! Tune_Access_And_Batching_Impl

   !! Reuse the transport decision, not MPI datatype or communicator ownership.
   subroutine Plan_Communication_Strategy(plan)
      type(Transfer_Plan), intent(inout) :: plan

      if (transport_policy >= 0) then
         !! Every explicitly selected transport also covers entry/export.
         !! Unsupported/count-limited paths retain their safe fallback.
         call Plan_Transport(plan, 0)
      else
         call Tune_Transfer(plan, planning_level)
      end if
   end subroutine ! Plan_Communication_Strategy

   subroutine Build_Transform_Plan( &
      plan, physical_layout, &
      dimensions, boundary_types, decomposition_type, internal_process_grid, &
      complex_input, generation_id)
      type(Transform_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: physical_layout
      integer, intent(in) :: dimensions, boundary_types(3), decomposition_type, internal_process_grid(3), generation_id
      logical, intent(in) :: complex_input

      call Build_Transform_Plan_Impl( &
         plan, physical_layout, &
         dimensions, boundary_types, decomposition_type, internal_process_grid, &
         complex_input, generation_id)
   end subroutine ! Build_Transform_Plan

   subroutine Build_Transform_Plan_Impl( &
      plan, physical_layout, &
      dimensions, boundary_types, decomposition_type, internal_process_grid, &
      complex_input, generation_id)
      type(Transform_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: physical_layout
      integer, intent(in) :: dimensions, boundary_types(3), decomposition_type, internal_process_grid(3), generation_id
      logical, intent(in) :: complex_input
      integer :: choice

      !! Both grades use the requested storage organization; the default preserves
      !! physical axis order. Tuning occurs within each required stage.
      choice = 3
      if (layout_policy >= 0) choice = layout_policy
      call Build_Candidate( &
         plan, physical_layout, &
         dimensions, boundary_types, decomposition_type, internal_process_grid, &
         complex_input, generation_id, choice)
   end subroutine ! Build_Transform_Plan_Impl

   !! Keep every forward layout transition paired with its reverse transition.
   !! plan%ready is set only after all active axes have prepared kernels.
   subroutine Build_Candidate( &
      plan, physical_layout, &
      dimensions, boundary_types, decomposition_type, internal_process_grid, &
      complex_input, generation_id, choice)
      type(Transform_Plan), intent(inout) :: plan
      type(Layout), intent(in) :: physical_layout
      integer, intent(in) :: dimensions, boundary_types(3), decomposition_type, internal_process_grid(3), &
         generation_id, choice
      logical, intent(in) :: complex_input
      type(Layout) :: current, lin, lout, owned, stage_in, stage_out
      integer :: axis, order(3), grid(3), n(3), rank, ierr, lines, np
      logical :: iscomplex

      call Release_Transform_Plan(plan)
      call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
      call MPI_Comm_size(ISFFT_Comm, np, ierr)
      plan%physical = physical_layout
      plan%dimensions = dimensions
      plan%generation = generation_id
      plan%layout_choice = choice
      current = physical_layout
      iscomplex = complex_input
      do axis = 1, dimensions
         order = [axis, 1, 2]
         select case (axis)
         case (1)
            order = [1, 2, 3]
         case (2)
            order = [2, 1, 3]
         case (3)
            order = [3, 1, 2]
         end select
         grid = internal_process_grid
         if (decomposition_type /= 0 .and. np > 1) then
            if (dimensions == 2) then
               grid = [1, np, 1]
               if (axis == 2) grid = [np, 1, 1]
            else if (decomposition_type == 1) then
               grid = [1, 1, np]
               if (axis == 3) grid = [1, np, 1]
            else
               select case (axis)
               case (1)
                  grid = internal_process_grid
               case (2)
                  grid = [internal_process_grid(2), 1, internal_process_grid(3)]
               case (3)
                  grid = [internal_process_grid(2), internal_process_grid(3), 1]
               end select
            end if
         end if
         lin = Make_Layout(current%global_shape, grid, order, iscomplex)
         if (decomposition_type == 0 .and. grid(axis) > 1) lin%replicated_axis = axis
         lout = lin
         plan%axes(axis)%transform_kind = boundary_types(axis) + 2
         if (boundary_types(axis) == 0 .and. .not. iscomplex) then
            plan%half_axis = axis
            plan%axes(axis)%transform_kind = 1
            lout%global_shape(axis) = current%global_shape(axis) / 2 + 1
            lout%complex_data = .true.
            iscomplex = .true.
         end if
         if (choice > 0 .and. plan%axes(axis)%transform_kind == 2) then
            lin%axis_order = [1, 2, 3]
            lout%axis_order = [1, 2, 3]
         end if
         stage_in = lin
         stage_out = lout
         if (choice > 0) then
            stage_in%axis_order = [1, 2, 3]
            stage_out%axis_order = [1, 2, 3]
         end if
         owned = stage_out
         if (decomposition_type == 0) then
            owned%replicated_axis = 0
            owned%axis_order = [1, 2, 3]
         end if
         plan%axes(axis)%physical_axis = axis
         plan%axes(axis)%transform_length = current%global_shape(axis)
         plan%axes(axis)%input = lin
         plan%axes(axis)%output = lout
         if (decomposition_type == 0 .and. grid(axis) > 1) then
            call Build_Block_Gather(plan%axes(axis)%before, current, stage_in)
            !! Execute_Block_Axis only needs before's gather communicator and the
            !! owned output geometry. Whole-field remaps below are never executed.
            plan%axes(axis)%after%source = stage_out
            plan%axes(axis)%after%destination = owned
         else
            !! Plan internal transposes and their reverse exchanges.
            call Build_Transfer(plan%axes(axis)%before, current, stage_in)
            call Build_Transfer(plan%axes(axis)%reverse_after, stage_in, current)
            call Build_Transfer(plan%axes(axis)%after, stage_out, owned)
            call Build_Transfer(plan%axes(axis)%reverse_before, owned, stage_out)
            call Build_Transfer(plan%axes(axis)%pack, stage_in, lin)
            call Build_Transfer(plan%axes(axis)%unpack, lout, stage_out)
            call Build_Transfer(plan%axes(axis)%reverse_pack, stage_out, lout)
            call Build_Transfer(plan%axes(axis)%reverse_unpack, lin, stage_in)
            call Plan_Transport(plan%axes(axis)%before, choice)
            call Plan_Transport(plan%axes(axis)%reverse_after, choice)
            call Plan_Transport(plan%axes(axis)%after, choice)
            call Plan_Transport(plan%axes(axis)%reverse_before, choice)
         end if
         n = Layout_Shape(lin, rank)
         if (decomposition_type == 0 .and. grid(axis) > 1) then
            !! Execution uses bounded tiles. Do not tune a full replicated field that is never executed.
            lines = int(product(int(n, ISFFT_INT64_KIND)) / plan%axes(axis)%transform_length)
            lines = Policy_Block_Capacity(plan%axes(axis)%transform_length, lines)
            plan%axes(axis)%tiled_block = .true.
            call Tune_Block_Width(plan%axes(axis))
         else
            lines = int(product(int(n, ISFFT_INT64_KIND)) / plan%axes(axis)%transform_length)
            call Build_Kernel_Plan(plan%axes(axis)%transform_length, plan%axes(axis)%transform_kind, lines, &
               plan%axes(axis)%kernel, &
               lin%complex_data)
            call Tune_Access_And_Batching(plan%axes(axis)%kernel, n, Axis_Position(lin%axis_order, axis))
         end if
         current = owned
      end do
      plan%spectral = current
      plan%ready = .true.
   end subroutine ! Build_Candidate

   !! Usage bits control prepared dependencies. Poisson multipliers belong to
   !! Poisson-containing scopes. Public FFT execution uses its transform plans.
   subroutine Initialize_Plans( &
      global_shape, prescribed_process_grid, internal_process_grid, &
      dimensions, boundary_types, decomposition_type, &
      domain_lengths, discretization_scheme, generation_id)
      integer, intent(in) :: global_shape(3), prescribed_process_grid(3), internal_process_grid(3), dimensions, &
         boundary_types(3), decomposition_type, discretization_scheme, generation_id
      real(kind=ISFFT_REAL_KIND), intent(in) :: domain_lengths(3)

      call Finalize_Plans()
      plan_dim = dimensions
      plan_bc = boundary_types
      plan_decomp = decomposition_type
      plan_generation = generation_id
      original_layout = Make_Layout(global_shape, prescribed_process_grid, [1, 2, 3], .false.)
      internal_layout = Make_Layout(global_shape, internal_process_grid, [1, 2, 3], .false.)
      if (btest(prepare_scope, 0)) then
         !! Plan redistribution into the internal decomposition and back.
         call Build_Transfer(redistribution_forward_plan, original_layout, internal_layout)
         call Plan_Communication_Strategy(redistribution_forward_plan)
         call Build_Transfer(redistribution_inverse_plan, internal_layout, original_layout)
         call Plan_Communication_Strategy(redistribution_inverse_plan)
      end if
      topology_ready = .true.
      if (btest(prepare_scope, 1)) then
         call Build_Transform_Plan(poisson_plan, internal_layout, dimensions, boundary_types, decomposition_type, &
            internal_process_grid, .false., generation_id)
         call Build_Poisson_Operator(poisson_symbol, poisson_plan%spectral, global_shape, boundary_types, dimensions, &
            domain_lengths, discretization_scheme)
      end if
   end subroutine ! Initialize_Plans

   subroutine Ensure_Public_Plan(is_r2c)
      logical, intent(in) :: is_r2c
      type(Layout) :: input
      integer :: which, mode, a
      type(Layout) :: wire

      which = 1
      if (is_r2c) which = 2
      if (public_plans(which)%ready) return
      input = internal_layout
      input%complex_data = .not. is_r2c
      if (is_r2c .and. poisson_plan%ready .and. all(plan_bc(:plan_dim) == 0)) then
         public_plans(which) = poisson_plan
         do a = 1, plan_dim
         end do
      else
         call Build_Transform_Plan(public_plans(which), input, plan_dim, [0, 0, 0], plan_decomp, &
            internal_layout%process_grid, &
            .not. is_r2c, plan_generation)
      end if
      do mode = 1, 2
         if (.not. btest(public_fft_modes(which), mode - 1)) cycle
         wire = internal_layout
         wire%complex_data = .true.
         wire%partition_shape = wire%global_shape
         if (is_r2c) wire%global_shape(1) = wire%global_shape(1) / 2 + 1
         if (mode == 1) wire = public_plans(which)%spectral
         public_layout(which, mode) = wire
         if (mode == 1) cycle
         call Build_Transfer(public_export(which, mode), public_plans(which)%spectral, wire)
         call Build_Transfer(public_import(which, mode), wire, public_plans(which)%spectral)
         call Plan_Communication_Strategy(public_export(which, mode))
         call Plan_Communication_Strategy(public_import(which, mode))
      end do
   end subroutine ! Ensure_Public_Plan

end module ISFFT_Planning
