module ISFFT_Repartitioning
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   implicit none

   public :: ISFFT_Repartitioning_Forward, ISFFT_Repartitioning_Inverse

   private

   integer, save :: original_npx0, original_npy0, original_npz0
   integer, save :: target_npx0, target_npy0, target_npz0
   integer, save :: original_i_offset(0:2048),original_j_offset(0:2048),original_k_offset(0:2048),&
      original_i_nn(0:2048),original_j_nn(0:2048),original_k_nn(0:2048)
   integer, save :: target_i_offset(0:2048),target_j_offset(0:2048),target_k_offset(0:2048),&
      target_i_nn(0:2048),target_j_nn(0:2048),target_k_nn(0:2048)

   logical, save :: If_need_repartitioning = .false.
   logical, save :: If_first_time_need_repartitioning = .true.
   logical, save :: If_repartitioning_initialized = .false.

contains

   subroutine need_repartitioning()

      if (.not.If_first_time_need_repartitioning) return

      if (np_size .eq. 1) then
         If_need_repartitioning = .false.
         If_first_time_need_repartitioning = .false.
         return
      endif

      select case (If_3dfft_decomp)
       case (0)  ! Block Decomposition
         If_need_repartitioning = .false.

       case (1)  ! Slab Decomposition
         if ((npx0 .ne. 1) .or. (npy0 .ne. 1) .or. (npz0 .ne. np_size)) then
            If_need_repartitioning = .true.
         end if

       case (2)  ! Pencil Decomposition
         if ((npx0 .ne. 1) .or. (npy .ne. 2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2)) .or. &
            (npz0 .ne. np_size/(2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2))) ) then
            If_need_repartitioning = .true.
         end if
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

      call Which_Target_Layout(target_npx0, target_npy0, target_npz0)

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
   subroutine Which_Target_Layout(npx0_target, npy0_target, npz0_target)
      integer, intent(out) :: npx0_target, npy0_target, npz0_target
      integer :: n1, n2, n3, ierr

      select case(If_3dfft_decomp)
       case (0)  ! Block Decomposition
         npx0_target = npx0
         npy0_target = npy0
         npz0_target = npz0

       case (1)  ! Slab Decomposition
         npx0_target = 1
         npy0_target = 1
         npz0_target = npx0*npy0*npz0

         n1 = nx_global/npx0_target
         n2 = ny_global/npy0_target
         n3 = nz_global/npz0_target

         if (np_size .gt. min(nx_global, ny_global, nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Slab Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

       case (2)  ! Pencil Decomposition
         npx0_target = 1
         npy0_target = 2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2)
         npz0_target = npx0*npy0*npz0/npy0_target

         n1 = nx_global/npx0_target
         n2 = ny_global/npy0_target
         n3 = nz_global/npz0_target

         if (np_size .gt. min(nx_global*ny_global, ny_global*nz_global, nx_global*nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Pencil Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

      end select

   end subroutine Which_Target_Layout
   !----------------------------------------------------------------------
   subroutine Repartitioning_Kernel(phi_inout, npx0_new, npy0_new, npz0_new)
      use ISFFT_MPI_Part

      real(kind=ISFFT_REAL_KIND), intent(inout), allocatable :: phi_inout(:,:,:)
      ! logical, intent(in) :: If_LAP
      integer, intent(in) :: npx0_new, npy0_new, npz0_new

      integer :: i, j, k, p, ierr, idx
      integer :: old_nx, old_ny, old_nz, old_npx, old_npy, old_npz
      integer :: old_npx0_local, old_npy0_local, old_npz0_local
      integer :: i_global, j_global, k_global
      integer :: i_local_new, j_local_new, k_local_new
      integer :: new_node_i, new_node_j, new_node_k
      integer :: target_id, proc_id, ix, iy, iz
      integer :: overlap_i_start, overlap_i_end
      integer :: overlap_j_start, overlap_j_end
      integer :: overlap_k_start, overlap_k_end
      logical :: is_to_target

      integer, allocatable :: old_i_offset(:), old_j_offset(:), old_k_offset(:)
      integer, allocatable :: old_i_nn(:), old_j_nn(:), old_k_nn(:)
      integer              :: sendcounts(0:np_size-1), recvcounts(0:np_size-1)
      integer              :: sdispls(0:np_size-1), rdispls(0:np_size-1)
      integer              :: pos(0:np_size-1)
      real(kind=ISFFT_REAL_KIND), allocatable :: sendbuf(:), recvbuf(:)


      is_to_target = (npx0_new == target_npx0 .and. npy0_new == target_npy0 .and. npz0_new == target_npz0)

      select case(is_to_target)

       case(.true.)

         allocate(old_i_offset(0:original_npx0-1), old_j_offset(0:original_npy0-1), old_k_offset(0:original_npz0-1))
         allocate(old_i_nn(0:original_npx0-1), old_j_nn(0:original_npy0-1), old_k_nn(0:original_npz0-1))

         old_npx0_local = original_npx0;
         old_npy0_local = original_npy0;
         old_npz0_local = original_npz0;

         old_npx = mod(my_id,old_npx0_local);
         old_npy = mod(my_id,old_npx0_local*old_npy0_local)/old_npx0_local;
         old_npz = my_id/(old_npx0_local*old_npy0_local);

         old_nx = nx_global/old_npx0_local;
         old_ny = ny_global/old_npy0_local;
         old_nz = nz_global/old_npz0_local;
         if(old_npx .lt. mod(nx_global,old_npx0_local)) old_nx=old_nx+1
         if(old_npy .lt. mod(ny_global,old_npy0_local)) old_ny=old_ny+1
         if(old_npz .lt. mod(nz_global,old_npz0_local)) old_nz=old_nz+1

         old_i_offset(0:original_npx0-1) = original_i_offset(0:original_npx0-1)
         old_j_offset(0:original_npy0-1) = original_j_offset(0:original_npy0-1)
         old_k_offset(0:original_npz0-1) = original_k_offset(0:original_npz0-1)
         old_i_nn(0:original_npx0-1) = original_i_nn(0:original_npx0-1)
         old_j_nn(0:original_npy0-1) = original_j_nn(0:original_npy0-1)
         old_k_nn(0:original_npz0-1) = original_k_nn(0:original_npz0-1)

       case(.false.)

         allocate(old_i_offset(0:target_npx0-1), old_j_offset(0:target_npy0-1), old_k_offset(0:target_npz0-1))
         allocate(old_i_nn(0:target_npx0-1), old_j_nn(0:target_npy0-1), old_k_nn(0:target_npz0-1))

         old_npx0_local = target_npx0;
         old_npy0_local = target_npy0;
         old_npz0_local = target_npz0;

         old_npx = mod(my_id,old_npx0_local);
         old_npy = mod(my_id,old_npx0_local*old_npy0_local)/old_npx0_local;
         old_npz = my_id/(old_npx0_local*old_npy0_local);

         old_nx = nx_global/old_npx0_local;
         old_ny = ny_global/old_npy0_local;
         old_nz = nz_global/old_npz0_local;
         if(old_npx .lt. mod(nx_global,old_npx0_local)) old_nx=old_nx+1
         if(old_npy .lt. mod(ny_global,old_npy0_local)) old_ny=old_ny+1
         if(old_npz .lt. mod(nz_global,old_npz0_local)) old_nz=old_nz+1

         old_i_offset(0:target_npx0-1) = target_i_offset(0:target_npx0-1)
         old_j_offset(0:target_npy0-1) = target_j_offset(0:target_npy0-1)
         old_k_offset(0:target_npz0-1) = target_k_offset(0:target_npz0-1)
         old_i_nn(0:target_npx0-1) = target_i_nn(0:target_npx0-1)
         old_j_nn(0:target_npy0-1) = target_j_nn(0:target_npy0-1)
         old_k_nn(0:target_npz0-1) = target_k_nn(0:target_npz0-1)

      end select

      call ISFFT_Part_change(npx0_new, npy0_new, npz0_new)

      sendcounts = 0

      do k = 1, old_nz
         k_global = old_k_offset(old_npz) + (k - 1)
         call ISFFT_get_k_node(k_global, new_node_k, k_local_new)
         do j = 1, old_ny
            j_global = old_j_offset(old_npy) + (j - 1)
            call ISFFT_get_j_node(j_global, new_node_j, j_local_new)
            do i = 1, old_nx
               i_global = old_i_offset(old_npx) + (i - 1)
               call ISFFT_get_i_node(i_global, new_node_i, i_local_new)
               target_id = ISFFT_get_id(new_node_i, new_node_j, new_node_k)
               sendcounts(target_id) = sendcounts(target_id) + 1
            end do
         end do
      end do

      call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      sdispls(0) = 0
      rdispls(0) = 0
      do p = 1, np_size-1
         sdispls(p) = sdispls(p-1) + sendcounts(p-1)
         rdispls(p) = rdispls(p-1) + recvcounts(p-1)
      end do

      allocate(sendbuf(0:sum(sendcounts)-1), source=0.0_ISFFT_REAL_KIND)
      allocate(recvbuf(0:sum(recvcounts)-1), source=0.0_ISFFT_REAL_KIND)

      pos(0:np_size-1) = sdispls(0:np_size-1)

      do k = 1, old_nz
         k_global = old_k_offset(old_npz) + (k - 1)
         call ISFFT_get_k_node(k_global, new_node_k, k_local_new)
         do j = 1, old_ny
            j_global = old_j_offset(old_npy) + (j - 1)
            call ISFFT_get_j_node(j_global, new_node_j, j_local_new)
            do i = 1, old_nx
               i_global = old_i_offset(old_npx) + (i - 1)
               call ISFFT_get_i_node(i_global, new_node_i, i_local_new)
               target_id = ISFFT_get_id(new_node_i, new_node_j, new_node_k)
               sendbuf(pos(target_id)) = phi_inout(i, j, k)
               pos(target_id) = pos(target_id) + 1
            end do
         end do
      end do

      call MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
         recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)


      if (allocated(phi_inout)) deallocate(phi_inout)
      allocate(phi_inout(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP), source=0.0_ISFFT_REAL_KIND)

      do proc_id = 0, np_size-1
         if (recvcounts(proc_id) .eq. 0) cycle

         iz = proc_id / (old_npx0_local * old_npy0_local)
         iy = (proc_id - iz * old_npx0_local * old_npy0_local) / old_npx0_local
         ix = mod(proc_id, old_npx0_local)

         overlap_i_start = max(old_i_offset(ix), i_offset(npx))
         overlap_i_end   = min(old_i_offset(ix) + old_i_nn(ix) - 1, i_offset(npx) + i_nn(npx) - 1)
         overlap_j_start = max(old_j_offset(iy), j_offset(npy))
         overlap_j_end   = min(old_j_offset(iy) + old_j_nn(iy) - 1, j_offset(npy) + j_nn(npy) - 1)
         overlap_k_start = max(old_k_offset(iz), k_offset(npz))
         overlap_k_end   = min(old_k_offset(iz) + old_k_nn(iz) - 1, k_offset(npz) + k_nn(npz) - 1)

         idx = rdispls(proc_id)
         do k_global = overlap_k_start, overlap_k_end
            do j_global = overlap_j_start, overlap_j_end
               do i_global = overlap_i_start, overlap_i_end
                  i_local_new = i_global - i_offset(npx) + 1
                  j_local_new = j_global - j_offset(npy) + 1
                  k_local_new = k_global - k_offset(npz) + 1
                  phi_inout(i_local_new, j_local_new, k_local_new) = recvbuf(idx)
                  idx = idx + 1
               end do
            end do
         end do
      end do

      deallocate(sendbuf, recvbuf, &
         old_i_offset, old_j_offset, old_k_offset, &
         old_i_nn, old_j_nn, old_k_nn)

   end subroutine Repartitioning_Kernel
   !----------------------------------------------------------------------
   subroutine Repartitioning_Interface(phi_inout, to_target)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: phi_inout(:,:,:)
      ! logical, intent(in) :: If_LAP
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


end module ISFFT_Repartitioning
