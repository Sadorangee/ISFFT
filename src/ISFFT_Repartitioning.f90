!-----------------------------------------------------------------------
! Repartitioning moves owned physical points only.  It is intentionally
! not coupled to halo/ghost-cell exchange because its primary use is to
! prepare data for complete FFT or other spectral-space operations.
! If valid halo values are required after repartitioning back to the
! caller's physical layout, update them explicitly through
! ISFFT_Halos_Exchange at the appropriate outer algorithm level.
!
! Allocatable inputs may use either of the two project conventions:
!   no halo : (1:nx,       1:ny,       1:nz)
!   halo    : (1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
! The output automatically preserves the detected input convention.
! LAP=0 is treated as the no-halo convention.
!-----------------------------------------------------------------------
module ISFFT_Repartitioning
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   implicit none

   public :: ISFFT_Repartitioning_Forward, ISFFT_Repartitioning_Inverse
   public :: ISFFT_Get_Target_Layout
   public :: ISFFT_Activate_Target_Layout, ISFFT_Activate_Original_Layout

   private

   integer, save :: original_npx0, original_npy0, original_npz0
   integer, save :: target_npx0, target_npy0, target_npz0
   integer, save :: original_i_offset(0:2048),original_j_offset(0:2048),original_k_offset(0:2048),&
      original_i_nn(0:2048),original_j_nn(0:2048),original_k_nn(0:2048)
   integer, save :: target_i_offset(0:2048),target_j_offset(0:2048),target_k_offset(0:2048),&
      target_i_nn(0:2048),target_j_nn(0:2048),target_k_nn(0:2048)

   type Repartitioning_Plan
      logical :: initialized = .false.
      integer :: send_total = 0, recv_total = 0
      integer, allocatable :: sendcounts(:), recvcounts(:)
      integer, allocatable :: sdispls(:), rdispls(:)
      integer, allocatable :: send_i_start(:), send_i_end(:)
      integer, allocatable :: send_j_start(:), send_j_end(:)
      integer, allocatable :: send_k_start(:), send_k_end(:)
      integer, allocatable :: recv_i_start(:), recv_i_end(:)
      integer, allocatable :: recv_j_start(:), recv_j_end(:)
      integer, allocatable :: recv_k_start(:), recv_k_end(:)
   end type Repartitioning_Plan

   logical, save :: If_need_repartitioning = .false.
   logical, save :: If_first_time_need_repartitioning = .true.
   logical, save :: If_repartitioning_initialized = .false.
   real(kind=ISFFT_REAL_KIND), allocatable, save :: repart_sendbuf(:), repart_recvbuf(:)
   type(Repartitioning_Plan), save, target :: repart_forward_plan, repart_inverse_plan

contains

   subroutine need_repartitioning()

      if (.not.If_first_time_need_repartitioning) return

      if (np_size .eq. 1) then
         If_need_repartitioning = .false.
         If_first_time_need_repartitioning = .false.
         return
      endif

      select case (Dimension)
      case (1)
         If_need_repartitioning = .false.
      case (2)
         select case (If_decomp)
         case (0)  ! 2D block keeps the external process grid.
            If_need_repartitioning = .false.
         case (1)  ! 2D slab target: x is local and y is distributed.
            If_need_repartitioning = (npx0 /= 1 .or. npy0 /= np_size .or. npz0 /= 1)
         case default
            call Repartitioning_Abort('Dimension=2 does not support pencil decomposition')
         end select
      case (3)
         select case (If_decomp)
         case (0)
            If_need_repartitioning = .false.
         case (1)
            If_need_repartitioning = (npx0 /= 1 .or. npy0 /= 1 .or. npz0 /= np_size)
         case (2)
            If_need_repartitioning = (npx0 /= 1 .or. &
               npy0 /= 2**((int(log(real(np_size))/log(2.0)) + 1) / 2) .or. &
               npz0 /= np_size/(2**((int(log(real(np_size))/log(2.0)) + 1) / 2)))
         end select
      end select

      If_first_time_need_repartitioning = .false.

   end subroutine need_repartitioning
   !----------------------------------------------------------------------
   subroutine Initialize_ISFFT_Repartitioning()
      integer :: k, ka, ierr

      if (If_repartitioning_initialized) then
         ! if (my_id == 0) then
         !    write(*,*) "Warning: Repartitioning is already initialized. Call Reset_ISFFT_Repartitioning first."
         ! endif
         return
      endif

      original_npx0 = npx0
      original_npy0 = npy0
      original_npz0 = npz0

      original_i_offset(0:original_npx0-1) = i_offset(0:original_npx0-1)
      original_j_offset(0:original_npy0-1) = j_offset(0:original_npy0-1)
      original_k_offset(0:original_npz0-1) = k_offset(0:original_npz0-1)
      original_i_nn(0:original_npx0-1) = i_nn(0:original_npx0-1)
      original_j_nn(0:original_npy0-1) = j_nn(0:original_npy0-1)
      original_k_nn(0:original_npz0-1) = k_nn(0:original_npz0-1)

      call ISFFT_Get_Target_Layout(target_npx0, target_npy0, target_npz0)

      do k=0,target_npx0-1
         ka=min(k,mod(nx_global,target_npx0))
         target_i_offset(k)=int(nx_global/target_npx0)*k+ka+1
         target_i_nn(k)=nx_global/target_npx0
         if(k .lt. mod(nx_global,target_npx0)) target_i_nn(k)=target_i_nn(k)+1
      enddo

      do k=0,target_npy0-1
         ka=min(k,mod(ny_global,target_npy0))
         target_j_offset(k)=int(ny_global/target_npy0)*k+ka+1
         target_j_nn(k)=ny_global/target_npy0
         if(k .lt. mod(ny_global,target_npy0)) target_j_nn(k)=target_j_nn(k)+1
      enddo

      do k=0,target_npz0-1
         ka=min(k,mod(nz_global,target_npz0))
         target_k_offset(k)=int(nz_global/target_npz0)*k+ka+1
         target_k_nn(k)=nz_global/target_npz0
         if(k .lt. mod(nz_global,target_npz0)) target_k_nn(k)=target_k_nn(k)+1
      enddo

      If_repartitioning_initialized = .true.

      call Build_Repartitioning_Plan(repart_forward_plan, &
         original_npx0, original_npy0, original_npz0, &
         original_i_offset, original_j_offset, original_k_offset, &
         original_i_nn, original_j_nn, original_k_nn, &
         target_npx0, target_npy0, target_npz0, &
         target_i_offset, target_j_offset, target_k_offset, &
         target_i_nn, target_j_nn, target_k_nn)

      call Build_Repartitioning_Plan(repart_inverse_plan, &
         target_npx0, target_npy0, target_npz0, &
         target_i_offset, target_j_offset, target_k_offset, &
         target_i_nn, target_j_nn, target_k_nn, &
         original_npx0, original_npy0, original_npz0, &
         original_i_offset, original_j_offset, original_k_offset, &
         original_i_nn, original_j_nn, original_k_nn)

      ! if (my_id == 0) then
      !    write(*,*) "----------------------------------------------------------------------"
      !    write(*,'(A, I0, A, I0, A, I0)') &
      !       ' Original partition (npx0 x npy0 x npz0): ', original_npx0, 'x', original_npy0, 'x', original_npz0
      !    write(*,'(A, I0, A, I0, A, I0)') &
      !       ' Target   partition (npx0 x npy0 x npz0): ', target_npx0, 'x', target_npy0, 'x', target_npz0
      !    write(*,*) "----------------------------------------------------------------------"
      ! endif

   end subroutine Initialize_ISFFT_Repartitioning
   !----------------------------------------------------------------------
   subroutine Build_Repartitioning_Plan(plan, old_npx0, old_npy0, old_npz0, &
      old_i_offset, old_j_offset, old_k_offset, old_i_nn, old_j_nn, old_k_nn, &
      new_npx0, new_npy0, new_npz0, new_i_offset, new_j_offset, new_k_offset, &
      new_i_nn, new_j_nn, new_k_nn)

      type(Repartitioning_Plan), intent(inout) :: plan
      integer, intent(in) :: old_npx0, old_npy0, old_npz0
      integer, intent(in) :: new_npx0, new_npy0, new_npz0
      integer, intent(in) :: old_i_offset(0:), old_j_offset(0:), old_k_offset(0:)
      integer, intent(in) :: old_i_nn(0:), old_j_nn(0:), old_k_nn(0:)
      integer, intent(in) :: new_i_offset(0:), new_j_offset(0:), new_k_offset(0:)
      integer, intent(in) :: new_i_nn(0:), new_j_nn(0:), new_k_nn(0:)

      integer :: p, ix, iy, iz
      integer :: old_npx_local, old_npy_local, old_npz_local
      integer :: new_npx_local, new_npy_local, new_npz_local
      integer :: old_i_start, old_i_end, old_j_start, old_j_end, old_k_start, old_k_end
      integer :: new_i_start, new_i_end, new_j_start, new_j_end, new_k_start, new_k_end

      if (allocated(plan%sendcounts)) then
         deallocate(plan%sendcounts, plan%recvcounts, plan%sdispls, plan%rdispls, &
            plan%send_i_start, plan%send_i_end, plan%send_j_start, plan%send_j_end, &
            plan%send_k_start, plan%send_k_end, plan%recv_i_start, plan%recv_i_end, &
            plan%recv_j_start, plan%recv_j_end, plan%recv_k_start, plan%recv_k_end)
      end if

      allocate(plan%sendcounts(0:np_size-1), plan%recvcounts(0:np_size-1))
      allocate(plan%sdispls(0:np_size-1), plan%rdispls(0:np_size-1))
      allocate(plan%send_i_start(0:np_size-1), plan%send_i_end(0:np_size-1))
      allocate(plan%send_j_start(0:np_size-1), plan%send_j_end(0:np_size-1))
      allocate(plan%send_k_start(0:np_size-1), plan%send_k_end(0:np_size-1))
      allocate(plan%recv_i_start(0:np_size-1), plan%recv_i_end(0:np_size-1))
      allocate(plan%recv_j_start(0:np_size-1), plan%recv_j_end(0:np_size-1))
      allocate(plan%recv_k_start(0:np_size-1), plan%recv_k_end(0:np_size-1))

      plan%sendcounts = 0
      plan%recvcounts = 0

      old_npx_local = mod(my_id, old_npx0)
      old_npy_local = mod(my_id, old_npx0 * old_npy0) / old_npx0
      old_npz_local = my_id / (old_npx0 * old_npy0)
      new_npx_local = mod(my_id, new_npx0)
      new_npy_local = mod(my_id, new_npx0 * new_npy0) / new_npx0
      new_npz_local = my_id / (new_npx0 * new_npy0)

      old_i_start = old_i_offset(old_npx_local)
      old_i_end = old_i_start + old_i_nn(old_npx_local) - 1
      old_j_start = old_j_offset(old_npy_local)
      old_j_end = old_j_start + old_j_nn(old_npy_local) - 1
      old_k_start = old_k_offset(old_npz_local)
      old_k_end = old_k_start + old_k_nn(old_npz_local) - 1

      new_i_start = new_i_offset(new_npx_local)
      new_i_end = new_i_start + new_i_nn(new_npx_local) - 1
      new_j_start = new_j_offset(new_npy_local)
      new_j_end = new_j_start + new_j_nn(new_npy_local) - 1
      new_k_start = new_k_offset(new_npz_local)
      new_k_end = new_k_start + new_k_nn(new_npz_local) - 1

      do p = 0, np_size-1
         iz = p / (new_npx0 * new_npy0)
         iy = (p - iz * new_npx0 * new_npy0) / new_npx0
         ix = mod(p, new_npx0)
         call Store_Overlap(old_i_start, old_i_end, old_j_start, old_j_end, old_k_start, old_k_end, &
            new_i_offset(ix), new_i_offset(ix) + new_i_nn(ix) - 1, &
            new_j_offset(iy), new_j_offset(iy) + new_j_nn(iy) - 1, &
            new_k_offset(iz), new_k_offset(iz) + new_k_nn(iz) - 1, &
            plan%send_i_start(p), plan%send_i_end(p), plan%send_j_start(p), plan%send_j_end(p), &
            plan%send_k_start(p), plan%send_k_end(p), plan%sendcounts(p))

         iz = p / (old_npx0 * old_npy0)
         iy = (p - iz * old_npx0 * old_npy0) / old_npx0
         ix = mod(p, old_npx0)
         call Store_Overlap(old_i_offset(ix), old_i_offset(ix) + old_i_nn(ix) - 1, &
            old_j_offset(iy), old_j_offset(iy) + old_j_nn(iy) - 1, &
            old_k_offset(iz), old_k_offset(iz) + old_k_nn(iz) - 1, &
            new_i_start, new_i_end, new_j_start, new_j_end, new_k_start, new_k_end, &
            plan%recv_i_start(p), plan%recv_i_end(p), plan%recv_j_start(p), plan%recv_j_end(p), &
            plan%recv_k_start(p), plan%recv_k_end(p), plan%recvcounts(p))
      end do

      call Prefix_Displacements(plan%sendcounts, plan%sdispls, plan%send_total)
      call Prefix_Displacements(plan%recvcounts, plan%rdispls, plan%recv_total)
      plan%initialized = .true.

   end subroutine Build_Repartitioning_Plan
   !----------------------------------------------------------------------
   subroutine Store_Overlap(a_i_start, a_i_end, a_j_start, a_j_end, a_k_start, a_k_end, &
      b_i_start, b_i_end, b_j_start, b_j_end, b_k_start, b_k_end, &
      out_i_start, out_i_end, out_j_start, out_j_end, out_k_start, out_k_end, count)

      integer, intent(in) :: a_i_start, a_i_end, a_j_start, a_j_end, a_k_start, a_k_end
      integer, intent(in) :: b_i_start, b_i_end, b_j_start, b_j_end, b_k_start, b_k_end
      integer, intent(out) :: out_i_start, out_i_end, out_j_start, out_j_end, out_k_start, out_k_end
      integer, intent(out) :: count

      out_i_start = max(a_i_start, b_i_start)
      out_i_end = min(a_i_end, b_i_end)
      out_j_start = max(a_j_start, b_j_start)
      out_j_end = min(a_j_end, b_j_end)
      out_k_start = max(a_k_start, b_k_start)
      out_k_end = min(a_k_end, b_k_end)

      if (out_i_start <= out_i_end .and. out_j_start <= out_j_end .and. out_k_start <= out_k_end) then
         count = (out_i_end - out_i_start + 1) * (out_j_end - out_j_start + 1) * &
            (out_k_end - out_k_start + 1)
      else
         count = 0
      end if

   end subroutine Store_Overlap
   !----------------------------------------------------------------------
   subroutine Prefix_Displacements(counts, displs, total)
      integer, intent(in) :: counts(0:)
      integer, intent(out) :: displs(0:)
      integer, intent(out) :: total
      integer :: p

      displs(0) = 0
      do p = 1, size(counts)-1
         displs(p) = displs(p-1) + counts(p-1)
      end do
      total = displs(size(counts)-1) + counts(size(counts)-1)

   end subroutine Prefix_Displacements
   !----------------------------------------------------------------------
   subroutine ISFFT_Get_Target_Layout(npx0_target, npy0_target, npz0_target)
      integer, intent(out) :: npx0_target, npy0_target, npz0_target
      integer :: n1, n2, n3, ierr

      select case (Dimension)
      case (1)
         npx0_target = 1
         npy0_target = 1
         npz0_target = 1
      case (2)
         select case (If_decomp)
         case (0)
            npx0_target = npx0
            npy0_target = npy0
            npz0_target = 1
         case (1)
            npx0_target = 1
            npy0_target = np_size
            npz0_target = 1
         case default
            call Repartitioning_Abort('Dimension=2 does not support pencil decomposition')
         end select
      case (3)
         select case(If_decomp)
         case (0)
            npx0_target = npx0
            npy0_target = npy0
            npz0_target = npz0
         case (1)
            npx0_target = 1
            npy0_target = 1
            npz0_target = np_size
            if (np_size > min(nx_global, ny_global, nz_global)) then
               call Repartitioning_Abort('The 3D slab decomposition is not suitable')
            end if
         case (2)
            npx0_target = 1
            npy0_target = 2**((int(log(real(np_size))/log(2.0)) + 1) / 2)
            npz0_target = np_size/npy0_target
            if (np_size > min(nx_global*ny_global, ny_global*nz_global, nx_global*nz_global)) then
               call Repartitioning_Abort('The 3D pencil decomposition is not suitable')
            end if
         end select
      end select

   end subroutine ISFFT_Get_Target_Layout

   subroutine Repartitioning_Abort(message)
      character(len=*), intent(in) :: message
      integer :: ierr

      if (my_id == 0) write(*,'(A,A)') 'ERROR: ', trim(message)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      stop 1
   end subroutine Repartitioning_Abort
   !----------------------------------------------------------------------
   subroutine Repartitioning_Kernel(phi_inout, npx0_new, npy0_new, npz0_new)
      use ISFFT_MPI_Part

      real(kind=ISFFT_REAL_KIND), intent(inout), allocatable :: phi_inout(:,:,:)
      integer, intent(in) :: npx0_new, npy0_new, npz0_new

      integer :: i, j, k, ierr, idx
      integer :: send_total, recv_total
      integer :: old_npx, old_npy, old_npz
      integer :: old_npx0_local, old_npy0_local, old_npz0_local
      integer :: old_nx_local, old_ny_local, old_nz_local
      integer :: j_global, k_global
      integer :: i_local_new, j_local_new, k_local_new
      integer :: target_id, proc_id
      integer :: overlap_i_start, overlap_i_end
      integer :: overlap_j_start, overlap_j_end
      integer :: overlap_k_start, overlap_k_end
      integer :: old_i_start, old_j_start, old_k_start
      logical :: is_to_target
      logical :: input_has_halo, input_has_no_halo
      type(Repartitioning_Plan), pointer :: plan

      integer              :: pos(0:np_size-1)


      is_to_target = (npx0_new == target_npx0 .and. npy0_new == target_npy0 .and. npz0_new == target_npz0)
      if (is_to_target) then
         plan => repart_forward_plan
      else
         plan => repart_inverse_plan
      end if

      select case(is_to_target)

       case(.true.)

         old_npx0_local = original_npx0;
         old_npy0_local = original_npy0;
         old_npz0_local = original_npz0;

         old_npx = mod(my_id,old_npx0_local);
         old_npy = mod(my_id,old_npx0_local*old_npy0_local)/old_npx0_local;
         old_npz = my_id/(old_npx0_local*old_npy0_local);

         old_i_start = original_i_offset(old_npx)
         old_j_start = original_j_offset(old_npy)
         old_k_start = original_k_offset(old_npz)
         old_nx_local = original_i_nn(old_npx)
         old_ny_local = original_j_nn(old_npy)
         old_nz_local = original_k_nn(old_npz)

       case(.false.)

         old_npx0_local = target_npx0;
         old_npy0_local = target_npy0;
         old_npz0_local = target_npz0;

         old_npx = mod(my_id,old_npx0_local);
         old_npy = mod(my_id,old_npx0_local*old_npy0_local)/old_npx0_local;
         old_npz = my_id/(old_npx0_local*old_npy0_local);

         old_i_start = target_i_offset(old_npx)
         old_j_start = target_j_offset(old_npy)
         old_k_start = target_k_offset(old_npz)
         old_nx_local = target_i_nn(old_npx)
         old_ny_local = target_j_nn(old_npy)
         old_nz_local = target_k_nn(old_npz)

      end select

      input_has_no_halo = &
         lbound(phi_inout,1) == 1 .and. ubound(phi_inout,1) == old_nx_local .and. &
         lbound(phi_inout,2) == 1 .and. ubound(phi_inout,2) == old_ny_local .and. &
         lbound(phi_inout,3) == 1 .and. ubound(phi_inout,3) == old_nz_local
      input_has_halo = LAP > 0 .and. &
         lbound(phi_inout,1) == 1-LAP .and. ubound(phi_inout,1) == old_nx_local+LAP .and. &
         lbound(phi_inout,2) == 1-LAP .and. ubound(phi_inout,2) == old_ny_local+LAP .and. &
         lbound(phi_inout,3) == 1-LAP .and. ubound(phi_inout,3) == old_nz_local+LAP

      if (.not. input_has_no_halo .and. .not. input_has_halo) then
         call Repartitioning_Abort('input bounds match neither the no-halo nor the LAP-halo convention')
      end if

      if (npx0 /= npx0_new .or. npy0 /= npy0_new .or. npz0 /= npz0_new) then
         call ISFFT_Part_change(npx0_new, npy0_new, npz0_new)
      end if

      send_total = plan%send_total
      recv_total = plan%recv_total
      if (.not. allocated(repart_sendbuf)) then
         allocate(repart_sendbuf(0:max(send_total, 1)-1))
      else if (size(repart_sendbuf) < send_total) then
         deallocate(repart_sendbuf)
         allocate(repart_sendbuf(0:send_total-1))
      end if
      if (.not. allocated(repart_recvbuf)) then
         allocate(repart_recvbuf(0:max(recv_total, 1)-1))
      else if (size(repart_recvbuf) < recv_total) then
         deallocate(repart_recvbuf)
         allocate(repart_recvbuf(0:recv_total-1))
      end if

      pos(0:np_size-1) = plan%sdispls(0:np_size-1)

      do target_id = 0, np_size-1
         if (plan%sendcounts(target_id) .eq. 0) cycle

         overlap_i_start = plan%send_i_start(target_id)
         overlap_i_end = plan%send_i_end(target_id)
         overlap_j_start = plan%send_j_start(target_id)
         overlap_j_end = plan%send_j_end(target_id)
         overlap_k_start = plan%send_k_start(target_id)
         overlap_k_end = plan%send_k_end(target_id)

         do k_global = overlap_k_start, overlap_k_end
            k = k_global - old_k_start + 1
            do j_global = overlap_j_start, overlap_j_end
               j = j_global - old_j_start + 1
               i = overlap_i_start - old_i_start + 1
               repart_sendbuf(pos(target_id):pos(target_id) + overlap_i_end - overlap_i_start) = &
                  phi_inout(i:i + overlap_i_end - overlap_i_start, j, k)
               pos(target_id) = pos(target_id) + overlap_i_end - overlap_i_start + 1
            end do
         end do
      end do

      call MPI_Alltoallv(repart_sendbuf, plan%sendcounts, plan%sdispls, ISFFT_DATA_TYPE, &
         repart_recvbuf, plan%recvcounts, plan%rdispls, ISFFT_DATA_TYPE, MPI_COMM_WORLD, ierr)


      if (allocated(phi_inout)) deallocate(phi_inout)
      if (input_has_halo) then
         allocate(phi_inout(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP), source=0.0_ISFFT_REAL_KIND)
      else
         allocate(phi_inout(1:nx, 1:ny, 1:nz), source=0.0_ISFFT_REAL_KIND)
      end if

      do proc_id = 0, np_size-1
         if (plan%recvcounts(proc_id) .eq. 0) cycle

         overlap_i_start = plan%recv_i_start(proc_id)
         overlap_i_end = plan%recv_i_end(proc_id)
         overlap_j_start = plan%recv_j_start(proc_id)
         overlap_j_end = plan%recv_j_end(proc_id)
         overlap_k_start = plan%recv_k_start(proc_id)
         overlap_k_end = plan%recv_k_end(proc_id)

         idx = plan%rdispls(proc_id)
         do k_global = overlap_k_start, overlap_k_end
            k_local_new = k_global - k_offset(npz) + 1
            do j_global = overlap_j_start, overlap_j_end
               i_local_new = overlap_i_start - i_offset(npx) + 1
               j_local_new = j_global - j_offset(npy) + 1
               phi_inout(i_local_new:i_local_new + overlap_i_end - overlap_i_start, &
                  j_local_new, k_local_new) = repart_recvbuf(idx:idx + overlap_i_end - overlap_i_start)
               idx = idx + overlap_i_end - overlap_i_start + 1
            end do
         end do
      end do

   end subroutine Repartitioning_Kernel
   !----------------------------------------------------------------------
   subroutine Repartitioning_Interface(phi_inout, to_target)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: phi_inout(:,:,:)
      logical, intent(in) :: to_target

      if (.not. If_repartitioning_initialized) then
         if (my_id == 0) then
            write(*,*) 'ERROR: Repartitioning is not initialized. Call Initialize_ISFFT_Repartitioning first.'
         endif
         return
      endif

      if (to_target) call Repartitioning_Kernel(phi_inout, target_npx0, target_npy0, target_npz0)

      if (.not. to_target) call Repartitioning_Kernel(phi_inout, original_npx0, original_npy0, original_npz0)

   end subroutine Repartitioning_Interface
   !----------------------------------------------------------------------




   !----------------------------------------------------------------------
   subroutine ISFFT_Repartitioning_Forward(phi_inout)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: phi_inout(:,:,:)

      call need_repartitioning()
      if (.not.If_need_repartitioning) return

      call Initialize_ISFFT_Repartitioning()

      call Repartitioning_Interface(phi_inout, .true.)

   end subroutine ISFFT_Repartitioning_Forward
   !----------------------------------------------------------------------
   subroutine ISFFT_Repartitioning_Inverse(phi_inout)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: phi_inout(:,:,:)

      call need_repartitioning()
      if (.not.If_need_repartitioning) return

      call Repartitioning_Interface(phi_inout, .false.)

   end subroutine ISFFT_Repartitioning_Inverse
   !----------------------------------------------------------------------
   subroutine ISFFT_Activate_Target_Layout()
      use ISFFT_MPI_Part, only: ISFFT_part_change

      call need_repartitioning()
      if (.not. If_need_repartitioning) return

      if (.not. If_repartitioning_initialized) call Initialize_ISFFT_Repartitioning()
      if (npx0 /= target_npx0 .or. npy0 /= target_npy0 .or. npz0 /= target_npz0) then
         call ISFFT_part_change(target_npx0, target_npy0, target_npz0)
      end if
   end subroutine ISFFT_Activate_Target_Layout
   !----------------------------------------------------------------------
   subroutine ISFFT_Activate_Original_Layout()
      use ISFFT_MPI_Part, only: ISFFT_part_change

      call need_repartitioning()
      if (.not. If_need_repartitioning) return

      if (.not. If_repartitioning_initialized) call Initialize_ISFFT_Repartitioning()
      if (npx0 /= original_npx0 .or. npy0 /= original_npy0 .or. npz0 /= original_npz0) then
         call ISFFT_part_change(original_npx0, original_npy0, original_npz0)
      end if
   end subroutine ISFFT_Activate_Original_Layout
   !----------------------------------------------------------------------


end module ISFFT_Repartitioning
