module ISFFT_LocalTranspose
    use ISFFT_Parameters
    use ISFFT_MPI_Constants
    use ISFFT_Buffer
    implicit none
 
 contains

 !----------------------------------------------------------------------
 subroutine local_transpose_12(phi_in, phi_im_in, phi_out, phi_im_out)
    use ISFFT_Buffer, only: opt_lt12_bi, opt_lt12_bj, opt_lt12_bk, opt_lt12_variant
    implicit none
    real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:),  phi_im_in(:,:,:)
    real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:), phi_im_out(:,:,:)

    integer :: n1, n2, n3
    integer :: ii, jj, kk, iend, jend, kend, i, j, k
    logical :: need_tune

    n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)
    ! need_tune = ((.not. If_lt12_tuned) .or. (n1 /= tuned_lt12_n1) .or. &
    !    (n2 /= tuned_lt12_n2) .or. (n3 /= tuned_lt12_n3))

    ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
    do kk = 1, n3, opt_lt12_bk
       kend = min(kk+opt_lt12_bk-1, n3)
       do jj = 1, n2, opt_lt12_bj
          jend = min(jj+opt_lt12_bj-1, n2)
          do ii = 1, n1, opt_lt12_bi
             iend = min(ii+opt_lt12_bi-1, n1)

             if (opt_lt12_variant == 0) then
                ! 0：k→i→j, sequential write
                do k = kk, kend
                   do i = ii, iend
                      !$omp simd
                      do j = jj, jend
                         phi_out(j,i,k)   = phi_in(i,j,k)
                         phi_im_out(j,i,k) = phi_im_in(i,j,k)
                      end do
                   end do
                end do
             else
                ! 1：k→j→i, sequential read
                do k = kk, kend
                   do j = jj, jend
                      !$omp simd
                      do i = ii, iend
                         phi_out(j,i,k)   = phi_in(i,j,k)
                         phi_im_out(j,i,k) = phi_im_in(i,j,k)
                      end do
                   end do
                end do
             end if

          end do
       end do
    end do
    ! $omp end parallel do

 end subroutine local_transpose_12
 !----------------------------------------------------------------------
 subroutine local_transpose_13(phi_in, phi_im_in, phi_out, phi_im_out)
    use ISFFT_Buffer, only: If_lt13_tuned, tuned_lt13_n1, tuned_lt13_n2, tuned_lt13_n3, &
       opt_lt13_bi, opt_lt13_bj, opt_lt13_bk, opt_lt13_variant
    use ISFFT_MPI_Constants, only: my_id
    implicit none
    real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:),  phi_im_in(:,:,:)
    real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:), phi_im_out(:,:,:)

    integer :: n1, n2, n3
    integer :: ii, jj, kk, iend, jend, kend, i, j, k
    logical :: need_tune

    n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)

    need_tune = ((.not. If_lt13_tuned) .or. (n1 /= tuned_lt13_n1) .or. &
       (n2 /= tuned_lt13_n2) .or. (n3 /= tuned_lt13_n3))

    ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
    do kk = 1, n3, opt_lt13_bk
       kend = min(kk+opt_lt13_bk-1, n3)
       do jj = 1, n2, opt_lt13_bj
          jend = min(jj+opt_lt13_bj-1, n2)
          do ii = 1, n1, opt_lt13_bi
             iend = min(ii+opt_lt13_bi-1, n1)

             if (opt_lt13_variant == 0) then
                ! 0：i→j→k, sequential write
                do i = ii, iend
                   do j = jj, jend
                      !$omp simd
                      do k = kk, kend
                         phi_out(k,j,i)   = phi_in(i,j,k)
                         phi_im_out(k,j,i) = phi_im_in(i,j,k)
                      end do
                   end do
                end do
             else
                ! 1：k→j→i, sequential read
                do k = kk, kend
                   do j = jj, jend
                      !$omp simd
                      do i = ii, iend
                         phi_out(k,j,i)   = phi_in(i,j,k)
                         phi_im_out(k,j,i) = phi_im_in(i,j,k)
                      end do
                   end do
                end do
             end if

          end do
       end do
    end do
    ! $omp end parallel do

 end subroutine local_transpose_13
 !----------------------------------------------------------------------


end module ISFFT_LocalTranspose