!----------------------------------------------------------------------
module SPLASH_Halos_Exchange
   use SPLASH_Parameters
   use SPLASH_MPI_Constants
   implicit none
   
   private
   public :: SPLASH_Update_Periodic_Boundary_xyz

contains

!----------------------------------------------------------------------
subroutine SPLASH_Update_Periodic_Boundary_xyz(f)
   implicit none
   real(kind=SPLASH_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
   call SPLASH_Update_Periodic_Boundary_x(f)
   call SPLASH_Update_Periodic_Boundary_y(f)
   call SPLASH_Update_Periodic_Boundary_z(f)
   return
end subroutine SPLASH_Update_Periodic_Boundary_xyz
 !----------------------------------------------------------------------
subroutine SPLASH_Update_Periodic_Boundary_x(f)
   use SPLASH_MPI_Part, only: SPLASH_get_id, SPLASH_mod
   implicit none
   integer :: Iperiodic1 = 1
   real(kind=SPLASH_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

   integer :: ierr, Status(MPI_STATUS_SIZE)
   integer :: hop, max_hops_left, max_hops_right
   integer :: rL, rR, proc_left, proc_right
   integer :: cnt_recv_left,  cnt_recv_right
   integer :: cnt_send_left,  cnt_send_right
   integer :: s, filled_left, filled_right
   integer :: i, j, k, k1
   integer :: tag_lr, tag_rl

   real(kind=SPLASH_REAL_KIND), allocatable :: send_right_buf(:), recv_left_buf(:)
   real(kind=SPLASH_REAL_KIND), allocatable :: send_left_buf(:),  recv_right_buf(:)

   if (LAP <= 0) return

   if (Iperiodic1 == 1) then
      max_hops_left  = npx0 - 1
      max_hops_right = npx0 - 1
   else
      max_hops_left  = npx
      max_hops_right = (npx0-1)-npx
   end if

   if (Iperiodic1 == 1 .and. npx0 == 1) then
      do k = 1, nz
         do j = 1, ny
            f(1-LAP:0, j, k) = f(nx-LAP+1:nx, j, k)
            f(nx+1:nx+LAP, j, k) = f(1:LAP, j, k)
         end do
      end do
      return
   end if

   ! if (Iperiodic1 /= 1 .and. npx0 == 1) needs special treatment

   allocate(send_right_buf(LAP*ny*nz), recv_left_buf(LAP*ny*nz))
   allocate(send_left_buf (LAP*ny*nz), recv_right_buf(LAP*ny*nz))

   filled_left  = 0
   filled_right = 0

   hop = 1
   do while ( (filled_left  < LAP .and. hop <= max_hops_left)  .or. &
      (filled_right < LAP .and. hop <= max_hops_right) )


      if (hop <= max_hops_left) then
         rL = SPLASH_mod(npx - hop, npx0)
         if (Iperiodic1 == 1 .or. npx - hop >= 0) then
            proc_left = SPLASH_get_id(rL, npy, npz)
         else
            proc_left = MPI_PROC_NULL
         end if
      else
         proc_left = MPI_PROC_NULL
      end if

      if (hop <= max_hops_right) then
         rR = SPLASH_mod(npx + hop, npx0)
         if (Iperiodic1 == 1 .or. npx + hop <= npx0-1) then
            proc_right = SPLASH_get_id(rR, npy, npz)
         else
            proc_right = MPI_PROC_NULL
         end if
      else
         proc_right = MPI_PROC_NULL
      end if


      cnt_recv_left  = 0
      if (proc_left  /= MPI_PROC_NULL .and. filled_left  < LAP) then
         cnt_recv_left = min( LAP - filled_left,  i_nn(rL) )
      end if

      cnt_recv_right = 0
      if (proc_right /= MPI_PROC_NULL .and. filled_right < LAP) then
         cnt_recv_right = min( LAP - filled_right, i_nn(rR) )
      end if


      cnt_send_right = 0
      if (proc_right /= MPI_PROC_NULL) then

         cnt_send_right = LAP
         do s=1,hop-1
            cnt_send_right = cnt_send_right - i_nn( SPLASH_mod(rR - s, npx0) )
         end do
         if (cnt_send_right < 0) cnt_send_right = 0
         cnt_send_right = min( cnt_send_right, i_nn(npx) )
      end if


      cnt_send_left = 0
      if (proc_left /= MPI_PROC_NULL) then

         cnt_send_left = LAP
         do s=1,hop-1
            cnt_send_left = cnt_send_left - i_nn( SPLASH_mod(rL + s, npx0) )
         end do
         if (cnt_send_left < 0) cnt_send_left = 0
         cnt_send_left = min( cnt_send_left, i_nn(npx) )
      end if


      tag_lr = 8000 + hop


      if (cnt_send_right > 0) then
         k1 = 0
         do k=1, nz
            do j=1, ny
               do i=nx-cnt_send_right+1, nx
                  k1 = k1 + 1
                  send_right_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_right_buf, cnt_send_right*ny*nz, SPLASH_DATA_TYPE, proc_right, tag_lr, &
         recv_left_buf,  cnt_recv_left *ny*nz, SPLASH_DATA_TYPE, proc_left,  tag_lr, &
         MPI_COMM_WORLD, Status, ierr )


      if (cnt_recv_left > 0) then
         k1 = 0
         do k=1, nz
            do j=1, ny
               do i=1, cnt_recv_left
                  k1 = k1 + 1

                  f( -filled_left - cnt_recv_left + i, j, k ) = recv_left_buf(k1)
               end do
            end do
         end do
         filled_left = filled_left + cnt_recv_left
      end if


      tag_rl = 9000 + hop


      if (cnt_send_left > 0) then
         k1 = 0
         do k=1, nz
            do j=1, ny
               do i=1, cnt_send_left
                  k1 = k1 + 1
                  send_left_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_left_buf,  cnt_send_left *ny*nz, SPLASH_DATA_TYPE, proc_left,  tag_rl, &
         recv_right_buf, cnt_recv_right*ny*nz, SPLASH_DATA_TYPE, proc_right, tag_rl, &
         MPI_COMM_WORLD, Status, ierr )


      if (cnt_recv_right > 0) then
         k1 = 0
         do k=1, nz
            do j=1, ny
               do i=1, cnt_recv_right
                  k1 = k1 + 1

                  f( nx + filled_right + i, j, k ) = recv_right_buf(k1)
               end do
            end do
         end do
         filled_right = filled_right + cnt_recv_right
      end if

      hop = hop + 1
   end do

   if (SPLASH_Barrier_level >= 1) call MPI_Barrier(MPI_COMM_WORLD, ierr)

   deallocate(send_right_buf, recv_left_buf, send_left_buf, recv_right_buf)
end subroutine SPLASH_Update_Periodic_Boundary_x
 !----------------------------------------------------------------------
subroutine SPLASH_Update_Periodic_Boundary_y(f)
   use SPLASH_MPI_Part, only: SPLASH_get_id, SPLASH_mod
   implicit none
   integer :: Iperiodic2 = 1
   real(kind=SPLASH_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

   integer :: ierr, Status(MPI_STATUS_SIZE)
   integer :: hop, max_hops_dn, max_hops_up
   integer :: rD, rU, proc_dn, proc_up
   integer :: cnt_recv_dn,  cnt_recv_up
   integer :: cnt_send_dn,  cnt_send_up
   integer :: s, filled_dn, filled_up
   integer :: i, j, k, k1
   integer :: tag_du, tag_ud

   real(kind=SPLASH_REAL_KIND), allocatable :: send_up_buf(:), recv_dn_buf(:)
   real(kind=SPLASH_REAL_KIND), allocatable :: send_dn_buf(:), recv_up_buf(:)

   if (LAP <= 0) return

   if (Iperiodic2 == 1) then
      max_hops_dn = npy0 - 1
      max_hops_up = npy0 - 1
   else
      max_hops_dn = npy
      max_hops_up = (npy0-1) - npy
   end if

   if (Iperiodic2 == 1 .and. npy0 == 1) then
      do k = 1, nz
         do i = 1, nx
            f(i, 1-LAP:0, k) = f(i, ny-LAP+1:ny, k)
            f(i, ny+1:ny+LAP, k) = f(i, 1:LAP, k)
         end do
      end do
      return
   end if

   ! if (Iperiodic2 /= 1 .and. npy0 == 1) needs special treatment

   allocate(send_up_buf(LAP*nx*nz), recv_dn_buf(LAP*nx*nz))
   allocate(send_dn_buf(LAP*nx*nz), recv_up_buf(LAP*nx*nz))

   filled_dn = 0
   filled_up = 0

   hop = 1
   do while ( (filled_dn < LAP .and. hop <= max_hops_dn) .or. &
      (filled_up < LAP .and. hop <= max_hops_up) )

      if (hop <= max_hops_dn) then
         rD = SPLASH_mod(npy - hop, npy0)
         if (Iperiodic2 == 1 .or. npy - hop >= 0) then
            proc_dn = SPLASH_get_id(npx, rD, npz)
         else
            proc_dn = MPI_PROC_NULL
         end if
      else
         proc_dn = MPI_PROC_NULL
      end if

      if (hop <= max_hops_up) then
         rU = SPLASH_mod(npy + hop, npy0)
         if (Iperiodic2 == 1 .or. npy + hop <= npy0-1) then
            proc_up = SPLASH_get_id(npx, rU, npz)
         else
            proc_up = MPI_PROC_NULL
         end if
      else
         proc_up = MPI_PROC_NULL
      end if

      cnt_recv_dn = 0
      if (proc_dn /= MPI_PROC_NULL .and. filled_dn < LAP) then
         cnt_recv_dn = min( LAP - filled_dn, j_nn(rD) )
      end if

      cnt_recv_up = 0
      if (proc_up /= MPI_PROC_NULL .and. filled_up < LAP) then
         cnt_recv_up = min( LAP - filled_up, j_nn(rU) )
      end if

      cnt_send_up = 0
      if (proc_up /= MPI_PROC_NULL) then
         cnt_send_up = LAP
         do s=1,hop-1
            cnt_send_up = cnt_send_up - j_nn( SPLASH_mod(rU - s, npy0) )
         end do
         if (cnt_send_up < 0) cnt_send_up = 0
         cnt_send_up = min( cnt_send_up, j_nn(npy) )
      end if

      cnt_send_dn = 0
      if (proc_dn /= MPI_PROC_NULL) then
         cnt_send_dn = LAP
         do s=1,hop-1
            cnt_send_dn = cnt_send_dn - j_nn( SPLASH_mod(rD + s, npy0) )
         end do
         if (cnt_send_dn < 0) cnt_send_dn = 0
         cnt_send_dn = min( cnt_send_dn, j_nn(npy) )
      end if

      tag_du = 8100 + hop
      if (cnt_send_up > 0) then
         k1 = 0
         do k=1, nz
            do j=ny-cnt_send_up+1, ny
               do i=1, nx
                  k1 = k1 + 1
                  send_up_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_up_buf,  cnt_send_up*nx*nz, SPLASH_DATA_TYPE, proc_up, tag_du, &
         recv_dn_buf,  cnt_recv_dn*nx*nz, SPLASH_DATA_TYPE, proc_dn, tag_du, &
         MPI_COMM_WORLD, Status, ierr )

      if (cnt_recv_dn > 0) then
         k1 = 0
         do k=1, nz
            do j=1, cnt_recv_dn
               do i=1, nx
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
         do k=1, nz
            do j=1, cnt_send_dn
               do i=1, nx
                  k1 = k1 + 1
                  send_dn_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_dn_buf,  cnt_send_dn*nx*nz, SPLASH_DATA_TYPE, proc_dn, tag_ud, &
         recv_up_buf,  cnt_recv_up*nx*nz, SPLASH_DATA_TYPE, proc_up, tag_ud, &
         MPI_COMM_WORLD, Status, ierr )

      if (cnt_recv_up > 0) then
         k1 = 0
         do k=1, nz
            do j=1, cnt_recv_up
               do i=1, nx
                  k1 = k1 + 1
                  f(i, ny + filled_up + j, k) = recv_up_buf(k1)
               end do
            end do
         end do
         filled_up = filled_up + cnt_recv_up
      end if

      hop = hop + 1
   end do

   if (SPLASH_Barrier_level >= 1) call MPI_Barrier(MPI_COMM_WORLD, ierr)

   deallocate(send_up_buf, recv_dn_buf, send_dn_buf, recv_up_buf)
end subroutine SPLASH_Update_Periodic_Boundary_y
 !----------------------------------------------------------------------
subroutine SPLASH_Update_Periodic_Boundary_z(f)
   use SPLASH_MPI_Part, only: SPLASH_get_id, SPLASH_mod
   implicit none
   integer :: Iperiodic3 = 1
   real(kind=SPLASH_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

   integer :: ierr, Status(MPI_STATUS_SIZE)
   integer :: hop, max_hops_bk, max_hops_fr
   integer :: rB, rF, proc_bk, proc_fr
   integer :: cnt_recv_bk,  cnt_recv_fr
   integer :: cnt_send_bk,  cnt_send_fr
   integer :: s, filled_bk, filled_fr
   integer :: i, j, k, k1
   integer :: tag_bf, tag_fb

   real(kind=SPLASH_REAL_KIND), allocatable :: send_fr_buf(:), recv_bk_buf(:)
   real(kind=SPLASH_REAL_KIND), allocatable :: send_bk_buf(:), recv_fr_buf(:)

   if (LAP <= 0) return

   if (Iperiodic3 == 1) then
      max_hops_bk = npz0 - 1
      max_hops_fr = npz0 - 1
   else
      max_hops_bk = npz
      max_hops_fr = (npz0-1) - npz
   end if

   if (Iperiodic3 == 1 .and. npz0 == 1) then
      do j = 1, ny
         do i = 1, nx
            f(i, j, 1-LAP:0) = f(i, j, nz-LAP+1:nz)
            f(i, j, nz+1:nz+LAP) = f(i, j, 1:LAP)
         end do
      end do
      return
   end if

   ! if (Iperiodic3 /= 1 .and. npz0 == 1) needs special treatment

   allocate(send_fr_buf(LAP*nx*ny), recv_bk_buf(LAP*nx*ny))
   allocate(send_bk_buf(LAP*nx*ny), recv_fr_buf(LAP*nx*ny))

   filled_bk = 0
   filled_fr = 0

   hop = 1
   do while ( (filled_bk < LAP .and. hop <= max_hops_bk) .or. &
      (filled_fr < LAP .and. hop <= max_hops_fr) )

      if (hop <= max_hops_bk) then
         rB = SPLASH_mod(npz - hop, npz0)
         if (Iperiodic3 == 1 .or. npz - hop >= 0) then
            proc_bk = SPLASH_get_id(npx, npy, rB)
         else
            proc_bk = MPI_PROC_NULL
         end if
      else
         proc_bk = MPI_PROC_NULL
      end if

      if (hop <= max_hops_fr) then
         rF = SPLASH_mod(npz + hop, npz0)
         if (Iperiodic3 == 1 .or. npz + hop <= npz0-1) then
            proc_fr = SPLASH_get_id(npx, npy, rF)
         else
            proc_fr = MPI_PROC_NULL
         end if
      else
         proc_fr = MPI_PROC_NULL
      end if

      cnt_recv_bk = 0
      if (proc_bk /= MPI_PROC_NULL .and. filled_bk < LAP) then
         cnt_recv_bk = min( LAP - filled_bk, k_nn(rB) )
      end if

      cnt_recv_fr = 0
      if (proc_fr /= MPI_PROC_NULL .and. filled_fr < LAP) then
         cnt_recv_fr = min( LAP - filled_fr, k_nn(rF) )
      end if

      cnt_send_fr = 0
      if (proc_fr /= MPI_PROC_NULL) then
         cnt_send_fr = LAP
         do s=1,hop-1
            cnt_send_fr = cnt_send_fr - k_nn( SPLASH_mod(rF - s, npz0) )
         end do
         if (cnt_send_fr < 0) cnt_send_fr = 0
         cnt_send_fr = min( cnt_send_fr, k_nn(npz) )
      end if

      cnt_send_bk = 0
      if (proc_bk /= MPI_PROC_NULL) then
         cnt_send_bk = LAP
         do s=1,hop-1
            cnt_send_bk = cnt_send_bk - k_nn( SPLASH_mod(rB + s, npz0) )
         end do
         if (cnt_send_bk < 0) cnt_send_bk = 0
         cnt_send_bk = min( cnt_send_bk, k_nn(npz) )
      end if

      tag_bf = 8200 + hop
      if (cnt_send_fr > 0) then
         k1 = 0
         do k=nz-cnt_send_fr+1, nz
            do j=1, ny
               do i=1, nx
                  k1 = k1 + 1
                  send_fr_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_fr_buf,  cnt_send_fr*nx*ny, SPLASH_DATA_TYPE, proc_fr, tag_bf, &
         recv_bk_buf,  cnt_recv_bk*nx*ny, SPLASH_DATA_TYPE, proc_bk, tag_bf, &
         MPI_COMM_WORLD, Status, ierr )

      if (cnt_recv_bk > 0) then
         k1 = 0
         do k=1, cnt_recv_bk
            do j=1, ny
               do i=1, nx
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
         do k=1, cnt_send_bk
            do j=1, ny
               do i=1, nx
                  k1 = k1 + 1
                  send_bk_buf(k1) = f(i, j, k)
               end do
            end do
         end do
      end if

      call MPI_Sendrecv( send_bk_buf,  cnt_send_bk*nx*ny, SPLASH_DATA_TYPE, proc_bk, tag_fb, &
         recv_fr_buf,  cnt_recv_fr*nx*ny, SPLASH_DATA_TYPE, proc_fr, tag_fb, &
         MPI_COMM_WORLD, Status, ierr )

      if (cnt_recv_fr > 0) then
         k1 = 0
         do k=1, cnt_recv_fr
            do j=1, ny
               do i=1, nx
                  k1 = k1 + 1
                  f(i, j, nz + filled_fr + k) = recv_fr_buf(k1)
               end do
            end do
         end do
         filled_fr = filled_fr + cnt_recv_fr
      end if

      hop = hop + 1
   end do

   if (SPLASH_Barrier_level >= 1) call MPI_Barrier(MPI_COMM_WORLD, ierr)

   deallocate(send_fr_buf, recv_bk_buf, send_bk_buf, recv_fr_buf)
end subroutine SPLASH_Update_Periodic_Boundary_z
 !----------------------------------------------------------------------

end module SPLASH_Halos_Exchange
 !----------------------------------------------------------------------
