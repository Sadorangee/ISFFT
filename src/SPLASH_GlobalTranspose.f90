module SPLASH_GlobalTranspose
    use SPLASH_Parameters
    use SPLASH_MPI_Constants
    use SPLASH_Buffer
    implicit none
 
 contains

!----------------------------------------------------------------------
 subroutine global_transpose_12(phi_re_in, phi_im_in, phi_re_out, phi_im_out)
    real(kind=SPLASH_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
    real(kind=SPLASH_REAL_KIND), intent(out) :: phi_re_out(:,:,:), phi_im_out(:,:,:)

    integer :: i, j, k, l, p, ierr, proc_id
    integer :: sendcounts(0:np_size-1), recvcounts(0:np_size-1)
    integer :: sdispls(0:np_size-1), rdispls(0:np_size-1)
    integer :: n1, n2, n3, idx, n2_over
    integer :: x_proc, y_proc, z_proc, proc_x, proc_y, proc_z
    integer :: send_size, recv_size, send_count_per_proc

    n1 = size(phi_re_in, 1)
    n2 = size(phi_re_in, 2)
    n3 = size(phi_re_in, 3)

    x_proc = mod(my_id, npx0)
    y_proc = mod(my_id, npx0*npy0)/npx0
    z_proc = my_id / (npx0 * npy0)

    sendcounts = 0
    recvcounts = 0
    sdispls = 0
    rdispls = 0

    send_count_per_proc = n1 * n2 * n3 / npy0

    do p = 0, np_size-1
       proc_x = mod(p, npx0)
       proc_y = mod(p,npx0*npy0)/npx0
       proc_z = p / (npx0 * npy0)

       if (proc_z == z_proc .and. proc_x == x_proc) then
          sendcounts(p) = send_count_per_proc
          sdispls(p) = proc_y * send_count_per_proc
       else
          sendcounts(p) = 0
          sdispls(p) = 0
       end if
    end do

    ! allocate(sendbuf_3DFFT(sum(sendcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(sendbuf_im_3DFFT(sum(sendcounts)))

    idx = 1
    n2_over = n2/npy0
    do l = 1, npy0
       do k = 1, n3
          do j = 1, n2_over  !! Here requries the 2nd dimension of input array must be global size.
             do i = 1, n1
                sendbuf_3DFFT_complex(1,idx) = phi_re_in(i,j+(l-1)*n2_over,k)
                sendbuf_3DFFT_complex(2,idx) = phi_im_in(i,j+(l-1)*n2_over,k)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, &
       recvcounts, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ierr)

    rdispls(0) = 0
    do p = 1, np_size-1
       rdispls(p) = rdispls(p-1) + recvcounts(p-1)
    end do

    ! allocate(recvbuf_3DFFT(sum(recvcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(recvbuf_im_3DFFT(sum(recvcounts)))

    call MPI_Alltoallv(sendbuf_3DFFT_complex, sendcounts*2, sdispls*2, SPLASH_DATA_TYPE, &
       recvbuf_3DFFT_complex, recvcounts*2, rdispls*2, SPLASH_DATA_TYPE, &
       MPI_COMM_WORLD, ierr)

    idx = 1
    do l = 1, npy0
       do k = 1, n3
          do j = 1, n2_over
             do i = 1, n1

                phi_re_out(i+(l-1)*n1, j, k) = recvbuf_3DFFT_complex(1,idx)
                phi_im_out(i+(l-1)*n1, j, k) = recvbuf_3DFFT_complex(2,idx)
                idx = idx + 1
             end do
          end do
       end do
    end do

 end subroutine global_transpose_12
 !----------------------------------------------------------------------
 subroutine global_transpose_21(phi_re_in, phi_im_in, phi_re_out, phi_im_out)
    real(kind=SPLASH_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
    real(kind=SPLASH_REAL_KIND), intent(out) :: phi_re_out(:,:,:), phi_im_out(:,:,:)

    integer :: i, j, k, l, p, ierr, proc_id
    integer :: sendcounts(0:np_size-1), recvcounts(0:np_size-1)
    integer :: sdispls(0:np_size-1), rdispls(0:np_size-1)
    integer :: n1, n2, n3, idx, n1_over
    integer :: x_proc, y_proc, z_proc, proc_x, proc_y, proc_z
    integer :: send_size, recv_size, send_count_per_proc


    n1 = size(phi_re_in, 1)
    n2 = size(phi_re_in, 2)
    n3 = size(phi_re_in, 3)

    x_proc = mod(my_id, npx0)
    y_proc = mod(my_id, npx0*npy0)/npx0
    z_proc = my_id / (npx0 * npy0)

    sendcounts = 0
    recvcounts = 0
    sdispls = 0
    rdispls = 0

    send_count_per_proc = n1 * n2 * n3 / npy0

    do p = 0, np_size-1
       proc_x = mod(p, npx0)
       proc_y = mod(p,npx0*npy0)/npx0
       proc_z = p / (npx0 * npy0)

       if (proc_z == z_proc .and. proc_x == x_proc) then
          sendcounts(p) = send_count_per_proc
          sdispls(p) = proc_y * send_count_per_proc
       else
          sendcounts(p) = 0
          sdispls(p) = 0
       end if
    end do

    ! allocate(sendbuf_3DFFT(sum(sendcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(sendbuf_im_3DFFT(sum(sendcounts)))

    idx = 1
    n1_over = n1/npy0
    do l = 1, npy0
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1_over  !! Here requries the 1st dimension of input array must be global size.
                sendbuf_3DFFT_complex(1,idx) = phi_re_in(i+(l-1)*n1_over,j,k)
                sendbuf_3DFFT_complex(2,idx) = phi_im_in(i+(l-1)*n1_over,j,k)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, &
       recvcounts, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ierr)

    rdispls(0) = 0
    do p = 1, np_size-1
       rdispls(p) = rdispls(p-1) + recvcounts(p-1)
    end do

    ! allocate(recvbuf_3DFFT(sum(recvcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(recvbuf_im_3DFFT(sum(recvcounts)))

    call MPI_Alltoallv(sendbuf_3DFFT_complex, sendcounts*2, sdispls*2, SPLASH_DATA_TYPE, &
       recvbuf_3DFFT_complex, recvcounts*2, rdispls*2, SPLASH_DATA_TYPE, &
       MPI_COMM_WORLD, ierr)

    idx = 1
    do l = 1, npy0
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1_over
                phi_re_out(i, j+(l-1)*n2, k) = recvbuf_3DFFT_complex(1,idx)
                phi_im_out(i, j+(l-1)*n2, k) = recvbuf_3DFFT_complex(2,idx)
                idx = idx + 1
             end do
          end do
       end do
    end do

 end subroutine global_transpose_21
 !----------------------------------------------------------------------
 subroutine global_transpose_13(phi_re_in, phi_im_in, phi_re_out, phi_im_out)
    real(kind=SPLASH_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
    real(kind=SPLASH_REAL_KIND), intent(out) :: phi_re_out(:,:,:), phi_im_out(:,:,:)

    integer :: i, j, k, l, p, ierr, proc_id
    integer :: sendcounts(0:np_size-1), recvcounts(0:np_size-1)
    integer :: sdispls(0:np_size-1), rdispls(0:np_size-1)
    integer :: n1, n2, n3, idx, n3_over
    integer :: x_proc, y_proc, z_proc, proc_x, proc_y, proc_z
    integer :: send_size, recv_size, send_count_per_proc

    n1 = size(phi_re_in, 1)
    n2 = size(phi_re_in, 2)
    n3 = size(phi_re_in, 3)

    x_proc = mod(my_id, npx0)
    y_proc = mod(my_id, npx0*npy0)/npx0
    z_proc = my_id / (npx0 * npy0)

    sendcounts = 0
    recvcounts = 0
    sdispls = 0
    rdispls = 0

    send_count_per_proc = n1 * n2 * n3 / npz0

    do p = 0, np_size-1

       proc_x = mod(p, npx0)
       proc_y = mod(p, npx0*npy0)/npx0
       proc_z = p / (npx0 * npy0)

       if (proc_y == y_proc .and. proc_x == x_proc) then
          sendcounts(p) = send_count_per_proc
          sdispls(p) = proc_z * send_count_per_proc
       else
          sendcounts(p) = 0
          sdispls(p) = 0
       end if

    end do

    idx = 1
    n3_over = n3/npz0
    do l = 1, npz0
       do k = 1, n3_over      !! Here requries the 3rd dimension of input array must be global size.
          do j = 1, n2
             do i = 1, n1
                sendbuf_3DFFT_complex(1,idx) = phi_re_in(i, j, k+(l-1)*n3_over)
                sendbuf_3DFFT_complex(2,idx) = phi_im_in(i, j, k+(l-1)*n3_over)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, &
       recvcounts, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ierr)

    rdispls(0) = 0
    do p = 1, np_size-1
       rdispls(p) = rdispls(p-1) + recvcounts(p-1)
    end do

    call MPI_Alltoallv(sendbuf_3DFFT_complex, sendcounts*2, sdispls*2, SPLASH_DATA_TYPE, &
       recvbuf_3DFFT_complex, recvcounts*2, rdispls*2, SPLASH_DATA_TYPE, &
       MPI_COMM_WORLD, ierr)


    idx = 1
    do l = 1, npz0
       do k = 1, n3_over
          do j = 1, n2
             do i = 1, n1
                phi_re_out(i+(l-1)*n1, j, k) = recvbuf_3DFFT_complex(1,idx)
                phi_im_out(i+(l-1)*n1, j, k) = recvbuf_3DFFT_complex(2,idx)
                idx = idx + 1
             end do
          end do
       end do
    end do

 end subroutine global_transpose_13
 !----------------------------------------------------------------------
 subroutine global_transpose_31(phi_re_in, phi_im_in, phi_re_out, phi_im_out)
    real(kind=SPLASH_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
    real(kind=SPLASH_REAL_KIND), intent(out) :: phi_re_out(:,:,:), phi_im_out(:,:,:)

    integer :: i, j, k, l, p, ierr, proc_id
    integer :: sendcounts(0:np_size-1), recvcounts(0:np_size-1)
    integer :: sdispls(0:np_size-1), rdispls(0:np_size-1)
    integer :: n1, n2, n3, idx, n1_over
    integer :: x_proc, y_proc, z_proc, proc_x, proc_y, proc_z
    integer :: send_size, recv_size, send_count_per_proc

    n1 = size(phi_re_in, 1)
    n2 = size(phi_re_in, 2)
    n3 = size(phi_re_in, 3)

    x_proc = mod(my_id, npx0)
    y_proc = mod(my_id, npx0*npy0)/npx0
    z_proc = my_id / (npx0 * npy0)

    sendcounts = 0
    recvcounts = 0
    sdispls = 0
    rdispls = 0

    send_count_per_proc = n1 * n2 * n3 / npz0

    do p = 0, np_size-1

       proc_x = mod(p, npx0)
       proc_y = mod(p, npx0*npy0)/npx0
       proc_z = p / (npx0 * npy0)

       if (proc_y == y_proc .and. proc_x == x_proc) then
          sendcounts(p) = send_count_per_proc
          sdispls(p) = proc_z * send_count_per_proc
       else
          sendcounts(p) = 0
          sdispls(p) = 0
       end if

    end do

    ! allocate(sendbuf_3DFFT(sum(sendcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(sendbuf_im_3DFFT(sum(sendcounts)))

    idx = 1
    n1_over = n1/npz0
    do l = 1, npz0
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1_over  !! Here requries the 1st dimension of input array must be global size.
                sendbuf_3DFFT_complex(1,idx) = phi_re_in(i+(l-1)*n1_over, j, k)
                sendbuf_3DFFT_complex(2,idx) = phi_im_in(i+(l-1)*n1_over, j, k)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, &
       recvcounts, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ierr)

    rdispls(0) = 0
    do p = 1, np_size-1
       rdispls(p) = rdispls(p-1) + recvcounts(p-1)
    end do

    ! allocate(recvbuf_3DFFT(sum(recvcounts)))      !! If the size of buffer here is not nx*ny*nz, uncomment this line.
    ! allocate(recvbuf_im_3DFFT(sum(recvcounts)))

    call MPI_Alltoallv(sendbuf_3DFFT_complex, sendcounts*2, sdispls*2, SPLASH_DATA_TYPE, &
       recvbuf_3DFFT_complex, recvcounts*2, rdispls*2, SPLASH_DATA_TYPE, &
       MPI_COMM_WORLD, ierr)

    idx = 1
    do l = 1, npz0
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1_over
                phi_re_out(i, j, k+(l-1)*n3) = recvbuf_3DFFT_complex(1,idx)
                phi_im_out(i, j, k+(l-1)*n3) = recvbuf_3DFFT_complex(2,idx)
                idx = idx + 1
             end do
          end do
       end do
    end do

 end subroutine global_transpose_31
 !----------------------------------------------------------------------


end module SPLASH_GlobalTranspose