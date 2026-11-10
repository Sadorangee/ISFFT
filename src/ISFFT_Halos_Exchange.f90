!------------------------------------------------------------------------------
! ISFFT_Halos_Exchange
! Per-axis halo exchange and physical boundary filling.
!------------------------------------------------------------------------------
module ISFFT_Halos_Exchange
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_REAL_KIND, ISFFT_REAL_MPI_TYPE
   use ISFFT_Problem_State, only: boundary_types, halo_width, i_nn, i_offset, &
      ISFFT_Barrier_level, j_nn, j_offset, k_nn, &
      k_offset, npx, npx0, npy, &
      npy0, npz, npz0, nx, &
      nx_global, ny, ny_global, nz, &
      nz_global
   use ISFFT_Workspace, only: halo_work
   use ISFFT_MPI_Part, only: ISFFT_get_id, ISFFT_mod
   implicit none
   private
   public :: Fill_Physical_X_Halos, Fill_Physical_Y_Halos, Fill_Physical_Z_Halos, ISFFT_Update_Boundary_x, ISFFT_Update_Boundary_xy
   public :: ISFFT_Update_Boundary_xyz, ISFFT_Update_Boundary_y, ISFFT_Update_Boundary_z

contains

   !------------------------------------------------------------------------------
   !! dimensions-level boundary update entry points
   !------------------------------------------------------------------------------
   subroutine ISFFT_Update_Boundary_xyz(f)
      implicit none
      real(kind=ISFFT_REAL_KIND) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)

      call ISFFT_Update_Boundary_x(f, boundary_types(1))
      call ISFFT_Update_Boundary_y(f, boundary_types(2))
      call ISFFT_Update_Boundary_z(f, boundary_types(3))
      return
   end subroutine ISFFT_Update_Boundary_xyz

   subroutine ISFFT_Update_Boundary_xy(f)
      implicit none
      real(kind=ISFFT_REAL_KIND) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)

      call ISFFT_Update_Boundary_x(f, boundary_types(1))
      call ISFFT_Update_Boundary_y(f, boundary_types(2))
   end subroutine ISFFT_Update_Boundary_xy

   !------------------------------------------------------------------------------
   !! Per-axis halo exchange
   !------------------------------------------------------------------------------
   !! Keep x/y/z exchange tags, buffer roles and hop participation aligned.
   subroutine ISFFT_Update_Boundary_x(f, BC1)
      use ISFFT_MPI_Part, only: ISFFT_get_id, ISFFT_mod
      implicit none
      integer, intent(in) :: BC1
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)

      integer :: ierr, Status(MPI_STATUS_SIZE)
      integer :: communication_hops
      integer :: hop, max_hops_left, max_hops_right
      integer :: rL, rR, proc_left, proc_right
      integer :: cnt_recv_left, cnt_recv_right
      integer :: cnt_send_left, cnt_send_right
      integer :: s, filled_left, filled_right
      integer :: i, j, k, k1
      integer :: tag_lr, tag_rl

      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_right_buf(:), recv_left_buf(:)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_left_buf(:), recv_right_buf(:)

      if (halo_width <= 0) return

      if (BC1 == 0) then
         max_hops_left = npx0 - 1
         max_hops_right = npx0 - 1
      else
         max_hops_left = npx
         max_hops_right = (npx0 - 1) - npx
      end if

      if (BC1 == 0 .and. npx0 == 1) then
         do k = 1, nz
            do j = 1, ny
               f(1 - halo_width:0, j, k) = f(nx - halo_width + 1:nx, j, k)
               f(nx + 1:nx + halo_width, j, k) = f(1:halo_width, j, k)
            end do
         end do
         return
      else if (BC1 /= 0 .and. npx0 == 1) then
         call Fill_Physical_X_Halos(f, BC1)
         return
      end if

      if (.not. allocated(halo_work)) error stop 'halo workspace not prepared'
      if (size(halo_work, 1) < halo_width * ny * nz) error stop 'halo workspace too small'
      send_right_buf => halo_work(1:halo_width * ny * nz, 1)
      recv_left_buf => halo_work(1:halo_width * ny * nz, 2)
      send_left_buf => halo_work(1:halo_width * ny * nz, 3)
      recv_right_buf => halo_work(1:halo_width * ny * nz, 4)

      filled_left = 0
      filled_right = 0

      !! All ranks must keep participating while another rank may need their data.
      !! A rank that filled its own halos can still be a sender in a later hop.
      communication_hops = npx0 - 1
      if (minval(i_nn(0:npx0 - 1)) > 0) communication_hops = &
         min(communication_hops, (halo_width + minval(i_nn(0:npx0 - 1)) - 1) / minval(i_nn(0:npx0 - 1)))
      hop = 1
      do while (hop <= communication_hops)

         if (hop <= max_hops_left) then
            rL = ISFFT_mod(npx - hop, npx0)
            if (BC1 == 0 .or. npx - hop >= 0) then
               proc_left = ISFFT_get_id(rL, npy, npz)
            else
               proc_left = MPI_PROC_NULL
            end if
         else
            proc_left = MPI_PROC_NULL
         end if

         if (hop <= max_hops_right) then
            rR = ISFFT_mod(npx + hop, npx0)
            if (BC1 == 0 .or. npx + hop <= npx0 - 1) then
               proc_right = ISFFT_get_id(rR, npy, npz)
            else
               proc_right = MPI_PROC_NULL
            end if
         else
            proc_right = MPI_PROC_NULL
         end if

         cnt_recv_left = 0
         if (proc_left /= MPI_PROC_NULL .and. filled_left < halo_width) then
            cnt_recv_left = min(halo_width - filled_left, i_nn(rL))
         end if

         cnt_recv_right = 0
         if (proc_right /= MPI_PROC_NULL .and. filled_right < halo_width) then
            cnt_recv_right = min(halo_width - filled_right, i_nn(rR))
         end if

         cnt_send_right = 0
         if (proc_right /= MPI_PROC_NULL) then

            cnt_send_right = halo_width
            do s = 1, hop - 1
               cnt_send_right = cnt_send_right - i_nn(ISFFT_mod(rR - s, npx0))
            end do
            if (cnt_send_right < 0) cnt_send_right = 0
            cnt_send_right = min(cnt_send_right, i_nn(npx))
         end if

         cnt_send_left = 0
         if (proc_left /= MPI_PROC_NULL) then

            cnt_send_left = halo_width
            do s = 1, hop - 1
               cnt_send_left = cnt_send_left - i_nn(ISFFT_mod(rL + s, npx0))
            end do
            if (cnt_send_left < 0) cnt_send_left = 0
            cnt_send_left = min(cnt_send_left, i_nn(npx))
         end if

         tag_lr = 8000 + hop

         if (cnt_send_right > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, ny
                  do i = nx - cnt_send_right + 1, nx
                     k1 = k1 + 1
                     send_right_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_right_buf, cnt_send_right * ny * nz, ISFFT_REAL_MPI_TYPE, proc_right, tag_lr, &
            recv_left_buf, &
            cnt_recv_left * ny * nz, ISFFT_REAL_MPI_TYPE, proc_left, tag_lr, ISFFT_Comm, Status, ierr)

         if (cnt_recv_left > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, ny
                  do i = 1, cnt_recv_left
                     k1 = k1 + 1

                     f(-filled_left - cnt_recv_left + i, j, k) = recv_left_buf(k1)
                  end do
               end do
            end do
            filled_left = filled_left + cnt_recv_left
         end if

         tag_rl = 9000 + hop

         if (cnt_send_left > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, ny
                  do i = 1, cnt_send_left
                     k1 = k1 + 1
                     send_left_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_left_buf, cnt_send_left * ny * nz, ISFFT_REAL_MPI_TYPE, proc_left, tag_rl, &
            recv_right_buf, cnt_recv_right * ny * nz, ISFFT_REAL_MPI_TYPE, proc_right, tag_rl, ISFFT_Comm, Status, ierr)

         if (cnt_recv_right > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, ny
                  do i = 1, cnt_recv_right
                     k1 = k1 + 1

                     f(nx + filled_right + i, j, k) = recv_right_buf(k1)
                  end do
               end do
            end do
            filled_right = filled_right + cnt_recv_right
         end if

         hop = hop + 1
      end do

      if (ISFFT_Barrier_level >= 1) call MPI_Barrier(ISFFT_Comm, ierr)
      if (BC1 /= 0) call Fill_Physical_X_Halos(f, BC1)

   end subroutine ISFFT_Update_Boundary_x

   subroutine ISFFT_Update_Boundary_y(f, BC2)
      use ISFFT_MPI_Part, only: ISFFT_get_id, ISFFT_mod
      implicit none
      integer, intent(in) :: BC2
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)

      integer :: ierr, Status(MPI_STATUS_SIZE)
      integer :: communication_hops
      integer :: hop, max_hops_dn, max_hops_up
      integer :: rD, rU, proc_dn, proc_up
      integer :: cnt_recv_dn, cnt_recv_up
      integer :: cnt_send_dn, cnt_send_up
      integer :: s, filled_dn, filled_up
      integer :: i, j, k, k1
      integer :: tag_du, tag_ud

      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_up_buf(:), recv_dn_buf(:)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_dn_buf(:), recv_up_buf(:)

      if (halo_width <= 0) return

      if (BC2 == 0) then
         max_hops_dn = npy0 - 1
         max_hops_up = npy0 - 1
      else
         max_hops_dn = npy
         max_hops_up = (npy0 - 1) - npy
      end if

      if (BC2 == 0 .and. npy0 == 1) then
         do k = 1, nz
            do i = 1, nx
               f(i, 1 - halo_width:0, k) = f(i, ny - halo_width + 1:ny, k)
               f(i, ny + 1:ny + halo_width, k) = f(i, 1:halo_width, k)
            end do
         end do
         return
      else if (BC2 /= 0 .and. npy0 == 1) then
         call Fill_Physical_Y_Halos(f, BC2)
         return
      end if

      if (.not. allocated(halo_work)) error stop 'halo workspace not prepared'
      if (size(halo_work, 1) < halo_width * nx * nz) error stop 'halo workspace too small'
      send_up_buf => halo_work(1:halo_width * nx * nz, 1)
      recv_dn_buf => halo_work(1:halo_width * nx * nz, 2)
      send_dn_buf => halo_work(1:halo_width * nx * nz, 3)
      recv_up_buf => halo_work(1:halo_width * nx * nz, 4)

      filled_dn = 0
      filled_up = 0

      !! All ranks must keep participating while another rank may need their data.
      !! A rank that filled its own halos can still be a sender in a later hop.
      communication_hops = npy0 - 1
      if (minval(j_nn(0:npy0 - 1)) > 0) communication_hops = &
         min(communication_hops, (halo_width + minval(j_nn(0:npy0 - 1)) - 1) / minval(j_nn(0:npy0 - 1)))
      hop = 1
      do while (hop <= communication_hops)

         if (hop <= max_hops_dn) then
            rD = ISFFT_mod(npy - hop, npy0)
            if (BC2 == 0 .or. npy - hop >= 0) then
               proc_dn = ISFFT_get_id(npx, rD, npz)
            else
               proc_dn = MPI_PROC_NULL
            end if
         else
            proc_dn = MPI_PROC_NULL
         end if

         if (hop <= max_hops_up) then
            rU = ISFFT_mod(npy + hop, npy0)
            if (BC2 == 0 .or. npy + hop <= npy0 - 1) then
               proc_up = ISFFT_get_id(npx, rU, npz)
            else
               proc_up = MPI_PROC_NULL
            end if
         else
            proc_up = MPI_PROC_NULL
         end if

         cnt_recv_dn = 0
         if (proc_dn /= MPI_PROC_NULL .and. filled_dn < halo_width) then
            cnt_recv_dn = min(halo_width - filled_dn, j_nn(rD))
         end if

         cnt_recv_up = 0
         if (proc_up /= MPI_PROC_NULL .and. filled_up < halo_width) then
            cnt_recv_up = min(halo_width - filled_up, j_nn(rU))
         end if

         cnt_send_up = 0
         if (proc_up /= MPI_PROC_NULL) then
            cnt_send_up = halo_width
            do s = 1, hop - 1
               cnt_send_up = cnt_send_up - j_nn(ISFFT_mod(rU - s, npy0))
            end do
            if (cnt_send_up < 0) cnt_send_up = 0
            cnt_send_up = min(cnt_send_up, j_nn(npy))
         end if

         cnt_send_dn = 0
         if (proc_dn /= MPI_PROC_NULL) then
            cnt_send_dn = halo_width
            do s = 1, hop - 1
               cnt_send_dn = cnt_send_dn - j_nn(ISFFT_mod(rD + s, npy0))
            end do
            if (cnt_send_dn < 0) cnt_send_dn = 0
            cnt_send_dn = min(cnt_send_dn, j_nn(npy))
         end if

         tag_du = 8100 + hop
         if (cnt_send_up > 0) then
            k1 = 0
            do k = 1, nz
               do j = ny - cnt_send_up + 1, ny
                  do i = 1, nx
                     k1 = k1 + 1
                     send_up_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_up_buf, cnt_send_up * nx * nz, ISFFT_REAL_MPI_TYPE, proc_up, tag_du, recv_dn_buf, &
            cnt_recv_dn * nx * nz, ISFFT_REAL_MPI_TYPE, proc_dn, tag_du, ISFFT_Comm, Status, ierr)

         if (cnt_recv_dn > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, cnt_recv_dn
                  do i = 1, nx
                     k1 = k1 + 1
                     f(i, -filled_dn - cnt_recv_dn + j, k) = recv_dn_buf(k1)
                  end do
               end do
            end do
            filled_dn = filled_dn + cnt_recv_dn
         end if

         tag_ud = 9100 + hop
         if (cnt_send_dn > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, cnt_send_dn
                  do i = 1, nx
                     k1 = k1 + 1
                     send_dn_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_dn_buf, cnt_send_dn * nx * nz, ISFFT_REAL_MPI_TYPE, proc_dn, tag_ud, recv_up_buf, &
            cnt_recv_up * nx * nz, ISFFT_REAL_MPI_TYPE, proc_up, tag_ud, ISFFT_Comm, Status, ierr)

         if (cnt_recv_up > 0) then
            k1 = 0
            do k = 1, nz
               do j = 1, cnt_recv_up
                  do i = 1, nx
                     k1 = k1 + 1
                     f(i, ny + filled_up + j, k) = recv_up_buf(k1)
                  end do
               end do
            end do
            filled_up = filled_up + cnt_recv_up
         end if

         hop = hop + 1
      end do

      if (ISFFT_Barrier_level >= 1) call MPI_Barrier(ISFFT_Comm, ierr)
      if (BC2 /= 0) call Fill_Physical_Y_Halos(f, BC2)

   end subroutine ISFFT_Update_Boundary_y

   subroutine ISFFT_Update_Boundary_z(f, BC3)
      use ISFFT_MPI_Part, only: ISFFT_get_id, ISFFT_mod
      implicit none
      integer, intent(in) :: BC3
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)

      integer :: ierr, Status(MPI_STATUS_SIZE)
      integer :: communication_hops
      integer :: hop, max_hops_bk, max_hops_fr
      integer :: rB, rF, proc_bk, proc_fr
      integer :: cnt_recv_bk, cnt_recv_fr
      integer :: cnt_send_bk, cnt_send_fr
      integer :: s, filled_bk, filled_fr
      integer :: i, j, k, k1
      integer :: tag_bf, tag_fb

      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_fr_buf(:), recv_bk_buf(:)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: send_bk_buf(:), recv_fr_buf(:)

      if (halo_width <= 0) return

      if (BC3 == 0) then
         max_hops_bk = npz0 - 1
         max_hops_fr = npz0 - 1
      else
         max_hops_bk = npz
         max_hops_fr = (npz0 - 1) - npz
      end if

      if (BC3 == 0 .and. npz0 == 1) then
         do j = 1, ny
            do i = 1, nx
               f(i, j, 1 - halo_width:0) = f(i, j, nz - halo_width + 1:nz)
               f(i, j, nz + 1:nz + halo_width) = f(i, j, 1:halo_width)
            end do
         end do
         return
      else if (BC3 /= 0 .and. npz0 == 1) then
         call Fill_Physical_Z_Halos(f, BC3)
         return
      end if

      if (.not. allocated(halo_work)) error stop 'halo workspace not prepared'
      if (size(halo_work, 1) < halo_width * nx * ny) error stop 'halo workspace too small'
      send_fr_buf => halo_work(1:halo_width * nx * ny, 1)
      recv_bk_buf => halo_work(1:halo_width * nx * ny, 2)
      send_bk_buf => halo_work(1:halo_width * nx * ny, 3)
      recv_fr_buf => halo_work(1:halo_width * nx * ny, 4)

      filled_bk = 0
      filled_fr = 0

      !! All ranks must keep participating while another rank may need their data.
      !! A rank that filled its own halos can still be a sender in a later hop.
      communication_hops = npz0 - 1
      if (minval(k_nn(0:npz0 - 1)) > 0) communication_hops = &
         min(communication_hops, (halo_width + minval(k_nn(0:npz0 - 1)) - 1) / minval(k_nn(0:npz0 - 1)))
      hop = 1
      do while (hop <= communication_hops)

         if (hop <= max_hops_bk) then
            rB = ISFFT_mod(npz - hop, npz0)
            if (BC3 == 0 .or. npz - hop >= 0) then
               proc_bk = ISFFT_get_id(npx, npy, rB)
            else
               proc_bk = MPI_PROC_NULL
            end if
         else
            proc_bk = MPI_PROC_NULL
         end if

         if (hop <= max_hops_fr) then
            rF = ISFFT_mod(npz + hop, npz0)
            if (BC3 == 0 .or. npz + hop <= npz0 - 1) then
               proc_fr = ISFFT_get_id(npx, npy, rF)
            else
               proc_fr = MPI_PROC_NULL
            end if
         else
            proc_fr = MPI_PROC_NULL
         end if

         cnt_recv_bk = 0
         if (proc_bk /= MPI_PROC_NULL .and. filled_bk < halo_width) then
            cnt_recv_bk = min(halo_width - filled_bk, k_nn(rB))
         end if

         cnt_recv_fr = 0
         if (proc_fr /= MPI_PROC_NULL .and. filled_fr < halo_width) then
            cnt_recv_fr = min(halo_width - filled_fr, k_nn(rF))
         end if

         cnt_send_fr = 0
         if (proc_fr /= MPI_PROC_NULL) then
            cnt_send_fr = halo_width
            do s = 1, hop - 1
               cnt_send_fr = cnt_send_fr - k_nn(ISFFT_mod(rF - s, npz0))
            end do
            if (cnt_send_fr < 0) cnt_send_fr = 0
            cnt_send_fr = min(cnt_send_fr, k_nn(npz))
         end if

         cnt_send_bk = 0
         if (proc_bk /= MPI_PROC_NULL) then
            cnt_send_bk = halo_width
            do s = 1, hop - 1
               cnt_send_bk = cnt_send_bk - k_nn(ISFFT_mod(rB + s, npz0))
            end do
            if (cnt_send_bk < 0) cnt_send_bk = 0
            cnt_send_bk = min(cnt_send_bk, k_nn(npz))
         end if

         tag_bf = 8200 + hop
         if (cnt_send_fr > 0) then
            k1 = 0
            do k = nz - cnt_send_fr + 1, nz
               do j = 1, ny
                  do i = 1, nx
                     k1 = k1 + 1
                     send_fr_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_fr_buf, cnt_send_fr * nx * ny, ISFFT_REAL_MPI_TYPE, proc_fr, tag_bf, recv_bk_buf, &
            cnt_recv_bk * nx * ny, ISFFT_REAL_MPI_TYPE, proc_bk, tag_bf, ISFFT_Comm, Status, ierr)

         if (cnt_recv_bk > 0) then
            k1 = 0
            do k = 1, cnt_recv_bk
               do j = 1, ny
                  do i = 1, nx
                     k1 = k1 + 1
                     f(i, j, -filled_bk - cnt_recv_bk + k) = recv_bk_buf(k1)
                  end do
               end do
            end do
            filled_bk = filled_bk + cnt_recv_bk
         end if

         tag_fb = 9200 + hop
         if (cnt_send_bk > 0) then
            k1 = 0
            do k = 1, cnt_send_bk
               do j = 1, ny
                  do i = 1, nx
                     k1 = k1 + 1
                     send_bk_buf(k1) = f(i, j, k)
                  end do
               end do
            end do
         end if

         call MPI_Sendrecv(send_bk_buf, cnt_send_bk * nx * ny, ISFFT_REAL_MPI_TYPE, proc_bk, tag_fb, recv_fr_buf, &
            cnt_recv_fr * nx * ny, ISFFT_REAL_MPI_TYPE, proc_fr, tag_fb, ISFFT_Comm, Status, ierr)

         if (cnt_recv_fr > 0) then
            k1 = 0
            do k = 1, cnt_recv_fr
               do j = 1, ny
                  do i = 1, nx
                     k1 = k1 + 1
                     f(i, j, nz + filled_fr + k) = recv_fr_buf(k1)
                  end do
               end do
            end do
            filled_fr = filled_fr + cnt_recv_fr
         end if

         hop = hop + 1
      end do

      if (ISFFT_Barrier_level >= 1) call MPI_Barrier(ISFFT_Comm, ierr)
      if (BC3 /= 0) call Fill_Physical_Z_Halos(f, BC3)

   end subroutine ISFFT_Update_Boundary_z

   !------------------------------------------------------------------------------
   !! Exterior nonperiodic halo filling
   !------------------------------------------------------------------------------
   subroutine Fill_Physical_X_Halos(f, boundary_type)
      integer, intent(in) :: boundary_type
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)
      integer :: s, first_global

      !! USER NOTICE: exterior halos on nonperiodic axes are temporarily zero-filled.
      !! This is a TEMPORARY placeholder, NOT a rigorous Dirichlet/Neumann
      !! boundary closure. Users must replace it for their physical problem.
      !! Test global coordinates, including wide halos on non-endpoint ranks.
      !! Domain-owned boundary nodes and in-domain MPI halo values are preserved.
      if (boundary_type == 0) return
      first_global = i_offset(npx)
      do s = 1, halo_width
         if (first_global - s < 1) f(1 - s, 1:ny, 1:nz) = 0.0_ISFFT_REAL_KIND
         if (first_global + nx - 1 + s > nx_global) f(nx + s, 1:ny, 1:nz) = 0.0_ISFFT_REAL_KIND
      end do
   end subroutine Fill_Physical_X_Halos

   subroutine Fill_Physical_Y_Halos(f, boundary_type)
      integer, intent(in) :: boundary_type
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)
      integer :: s, first_global

      !! USER NOTICE: exterior halos on nonperiodic axes are temporarily zero-filled.
      !! This is a TEMPORARY placeholder, NOT a rigorous Dirichlet/Neumann
      !! boundary closure. Users must replace it for their physical problem.
      !! Test global coordinates, including wide halos on non-endpoint ranks.
      !! Domain-owned boundary nodes and in-domain MPI halo values are preserved.
      if (boundary_type == 0) return
      first_global = j_offset(npy)
      do s = 1, halo_width
         if (first_global - s < 1) f(1:nx, 1 - s, 1:nz) = 0.0_ISFFT_REAL_KIND
         if (first_global + ny - 1 + s > ny_global) f(1:nx, ny + s, 1:nz) = 0.0_ISFFT_REAL_KIND
      end do
   end subroutine Fill_Physical_Y_Halos

   subroutine Fill_Physical_Z_Halos(f, boundary_type)
      integer, intent(in) :: boundary_type
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)
      integer :: s, first_global

      !! USER NOTICE: exterior halos on nonperiodic axes are temporarily zero-filled.
      !! This is a TEMPORARY placeholder, NOT a rigorous Dirichlet/Neumann
      !! boundary closure. Users must replace it for their physical problem.
      !! Test global coordinates, including wide halos on non-endpoint ranks.
      !! Domain-owned boundary nodes and in-domain MPI halo values are preserved.
      if (boundary_type == 0) return
      first_global = k_offset(npz)
      do s = 1, halo_width
         if (first_global - s < 1) f(1:nx, 1:ny, 1 - s) = 0.0_ISFFT_REAL_KIND
         if (first_global + nz - 1 + s > nz_global) f(1:nx, 1:ny, nz + s) = 0.0_ISFFT_REAL_KIND
      end do
   end subroutine Fill_Physical_Z_Halos

end module ISFFT_Halos_Exchange
