module ISFFT_LocalTranspose
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   implicit none

contains

   !----------------------------------------------------------------------
   subroutine local_transpose_01_complex(phi_in, phi_im_in, phi_out, phi_im_out)
      use ISFFT_Buffer, only: opt_lt01_bi, opt_lt01_bj, opt_lt01_bk, opt_lt01_variant
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:),  phi_im_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:), phi_im_out(:,:,:)

      integer :: n1, n2, n3
      integer :: ii, jj, kk, iend, jend, kend, i, j, k
      logical :: need_tune

      n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)
      ! need_tune = ((.not. If_lt01_tuned) .or. (n1 /= tuned_lt01_n1) .or. &
      !    (n2 /= tuned_lt01_n2) .or. (n3 /= tuned_lt01_n3))

      ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
      do kk = 1, n3, opt_lt01_bk
         kend = min(kk+opt_lt01_bk-1, n3)
         do jj = 1, n2, opt_lt01_bj
            jend = min(jj+opt_lt01_bj-1, n2)
            do ii = 1, n1, opt_lt01_bi
               iend = min(ii+opt_lt01_bi-1, n1)

               if (opt_lt01_variant == 0) then
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

   end subroutine local_transpose_01_complex
   !----------------------------------------------------------------------
   subroutine local_transpose_02_complex(phi_in, phi_im_in, phi_out, phi_im_out)
      use ISFFT_Buffer, only: If_lt02_tuned, tuned_lt02_n1, tuned_lt02_n2, tuned_lt02_n3, &
         opt_lt02_bi, opt_lt02_bj, opt_lt02_bk, opt_lt02_variant
      use ISFFT_MPI_Constants, only: my_id
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:),  phi_im_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:), phi_im_out(:,:,:)

      integer :: n1, n2, n3
      integer :: ii, jj, kk, iend, jend, kend, i, j, k
      logical :: need_tune

      n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)

      need_tune = ((.not. If_lt02_tuned) .or. (n1 /= tuned_lt02_n1) .or. &
         (n2 /= tuned_lt02_n2) .or. (n3 /= tuned_lt02_n3))

      ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
      do kk = 1, n3, opt_lt02_bk
         kend = min(kk+opt_lt02_bk-1, n3)
         do jj = 1, n2, opt_lt02_bj
            jend = min(jj+opt_lt02_bj-1, n2)
            do ii = 1, n1, opt_lt02_bi
               iend = min(ii+opt_lt02_bi-1, n1)

               if (opt_lt02_variant == 0) then
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

   end subroutine local_transpose_02_complex
   !----------------------------------------------------------------------
   subroutine local_transpose_01_real(phi_in, phi_out)
      use ISFFT_Buffer, only: opt_lt01_bi, opt_lt01_bj, opt_lt01_bk, opt_lt01_variant
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:)

      integer :: n1, n2, n3
      integer :: ii, jj, kk, iend, jend, kend, i, j, k
      logical :: need_tune

      n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)
      ! need_tune = ((.not. If_lt01_tuned) .or. (n1 /= tuned_lt01_n1) .or. &
      !    (n2 /= tuned_lt01_n2) .or. (n3 /= tuned_lt01_n3))

      ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
      do kk = 1, n3, opt_lt01_bk
         kend = min(kk+opt_lt01_bk-1, n3)
         do jj = 1, n2, opt_lt01_bj
            jend = min(jj+opt_lt01_bj-1, n2)
            do ii = 1, n1, opt_lt01_bi
               iend = min(ii+opt_lt01_bi-1, n1)

               if (opt_lt01_variant == 0) then
                  ! 0：k→i→j, sequential write
                  do k = kk, kend
                     do i = ii, iend
                        !$omp simd
                        do j = jj, jend
                           phi_out(j,i,k) = phi_in(i,j,k)
                        end do
                     end do
                  end do
               else
                  ! 1：k→j→i, sequential read
                  do k = kk, kend
                     do j = jj, jend
                        !$omp simd
                        do i = ii, iend
                           phi_out(j,i,k) = phi_in(i,j,k)
                        end do
                     end do
                  end do
               end if

            end do
         end do
      end do
      ! $omp end parallel do

   end subroutine local_transpose_01_real
   !----------------------------------------------------------------------
   subroutine local_transpose_02_real(phi_in, phi_out)
      use ISFFT_Buffer, only: If_lt02_tuned, tuned_lt02_n1, tuned_lt02_n2, tuned_lt02_n3, &
         opt_lt02_bi, opt_lt02_bj, opt_lt02_bk, opt_lt02_variant
      use ISFFT_MPI_Constants, only: my_id
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(out)   :: phi_out(:,:,:)

      integer :: n1, n2, n3
      integer :: ii, jj, kk, iend, jend, kend, i, j, k
      logical :: need_tune

      n1 = size(phi_in,1); n2 = size(phi_in,2); n3 = size(phi_in,3)

      need_tune = ((.not. If_lt02_tuned) .or. (n1 /= tuned_lt02_n1) .or. &
         (n2 /= tuned_lt02_n2) .or. (n3 /= tuned_lt02_n3))

      ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
      do kk = 1, n3, opt_lt02_bk
         kend = min(kk+opt_lt02_bk-1, n3)
         do jj = 1, n2, opt_lt02_bj
            jend = min(jj+opt_lt02_bj-1, n2)
            do ii = 1, n1, opt_lt02_bi
               iend = min(ii+opt_lt02_bi-1, n1)

               if (opt_lt02_variant == 0) then
                  ! 0：i→j→k, sequential write
                  do i = ii, iend
                     do j = jj, jend
                        !$omp simd
                        do k = kk, kend
                           phi_out(k,j,i) = phi_in(i,j,k)
                        end do
                     end do
                  end do
               else
                  ! 1：k→j→i, sequential read
                  do k = kk, kend
                     do j = jj, jend
                        !$omp simd
                        do i = ii, iend
                           phi_out(k,j,i) = phi_in(i,j,k)
                        end do
                     end do
                  end do
               end if

            end do
         end do
      end do
      ! $omp end parallel do

   end subroutine local_transpose_02_real
   !----------------------------------------------------------------------



end module ISFFT_LocalTranspose
