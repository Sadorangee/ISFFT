module ISFFT_LT_Tuning
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   implicit none

contains

   !----------------------------------------------------------------------
   subroutine Tune_LocalTranspose(transpose_type)
      use, intrinsic :: iso_fortran_env, only: int64, real64
      implicit none
      integer, intent(in) :: transpose_type  ! 12 / 13

      integer :: n1, n2, n3, ierr
      real(kind=ISFFT_REAL_KIND), allocatable :: test_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), allocatable :: test_out(:,:,:)

      ! Target grid size
      if(If_3dfft_decomp .eq. 1) then
         n1 = nx_global/1
         n2 = ny_global/1
         n3 = nz_global/(np_size)

         if (np_size .gt. min(nx_global, ny_global, nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Slab Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

      else if(If_3dfft_decomp .eq. 2) then
         n1 = nx_global/1
         n2 = ny_global/(2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2))
         n3 = nz_global/(np_size/(2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2)))

         if (np_size .gt. min(nx_global*ny_global, ny_global*nz_global, nx_global*nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Pencil Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

      end if

      ! Allocate test arrays based on transpose type
      if (transpose_type == 12) then
         ! 12 transpose: exchange 1st and 2nd dimensions
         allocate(test_in(n1,n2,n3))
         allocate(test_out(n2,n1,n3))
      else if (transpose_type == 13) then
         ! 13 transpose: exchange 1st and 3rd dimensions
         allocate(test_in(n1,n2,n3))
         allocate(test_out(n3,n2,n1))
      else
         if (my_id == 0) write(*,*) "Error: Invalid transpose_type = ", transpose_type
         return
      end if

      ! Initialize test data
      call random_number(test_in)

      ! Perform tuning
      call Tune_Transpose_BS(test_in, test_out, &
         n1, n2, n3, transpose_type)

      ! Update tuning status
      if (transpose_type == 12) then
         If_lt12_tuned = .true.
         tuned_lt12_n1 = n1
         tuned_lt12_n2 = n2
         tuned_lt12_n3 = n3
         if (my_id == 0) then
            write(*,*) "Local Transpose 12 tuning completed:"
            write(*,*) "  Block sizes: bi=", opt_lt12_bi, " bj=", opt_lt12_bj, " bk=", opt_lt12_bk
            write(*,*) "  Variant: ", opt_lt12_variant
         end if
      else if (transpose_type == 13) then
         If_lt13_tuned = .true.
         tuned_lt13_n1 = n1
         tuned_lt13_n2 = n2
         tuned_lt13_n3 = n3
         if (my_id == 0) then
            write(*,*) "Local Transpose 13 tuning completed:"
            write(*,*) "  Block sizes: bi=", opt_lt13_bi, " bj=", opt_lt13_bj, " bk=", opt_lt13_bk
            write(*,*) "  Variant: ", opt_lt13_variant
         end if
      end if

      ! Deallocate test arrays
      deallocate(test_in, test_out)

   contains
      !------------------------------------------------------------
      subroutine Tune_Transpose_BS(ain, aout, dim1,dim2,dim3, trans_type)
         real(kind=ISFFT_REAL_KIND), intent(in)  :: ain(:,:,:)
         real(kind=ISFFT_REAL_KIND), intent(out) :: aout(:,:,:)
         integer, intent(in)  :: dim1,dim2,dim3, trans_type
         integer :: bi, bj, bk, variant

         ! Set candidate parameters for specific environments
         integer, parameter :: n_bi_candidates_12(*) = [32,48,64,96,128]
         integer, parameter :: n_bj_candidates_12(*) = [4,16,32]
         integer, parameter :: n_bk_candidates_12(*) = [1,2,4]

         integer, parameter :: n_bi_candidates_13(*) = [4,8,24,48]
         integer, parameter :: n_bj_candidates_13(*) = [1,2,4,8]
         integer, parameter :: n_bk_candidates_13(*) = [4,16,32]

         integer, parameter :: variants(*) = [0,1]

         integer :: ci, cj, ck, cv
         integer :: c_bi, c_bj, c_bk, c_var
         integer :: s1, s2, s3
         real(real64) :: best_t, t

         best_t = huge(1.0_real64)

         ! Set default values
         if (trans_type == 12) then
            bi=32; bj=32; bk=4; variant=0
            s1 = min(dim1, 128)
            s2 = min(dim2, 128)
            s3 = min(dim3, 64)
         else
            bi=32; bj=4; bk=32; variant=0
            s1 = min(dim1, 128)
            s2 = min(dim2, 64)
            s3 = min(dim3, 128)
         end if

         do cv = 1, size(variants)
            c_var = variants(cv)

            if (trans_type == 12) then
               do ci = 1, size(n_bi_candidates_12)
                  c_bi = n_bi_candidates_12(ci)
                  if (c_bi > s1) cycle
                  do cj = 1, size(n_bj_candidates_12)
                     c_bj = n_bj_candidates_12(cj)
                     if (c_bj > s2) cycle
                     do ck = 1, size(n_bk_candidates_12)
                        c_bk = n_bk_candidates_12(ck)
                        if (c_bk > s3) cycle

                        ! L1 cache limit
                        if ( c_bi*c_bj*c_bk*ISFFT_ELEMSIZE() > 24576*3 ) cycle

                        t = Bench_Transpose_Once(ain, aout, s1,s2,s3, c_bi,c_bj,c_bk, c_var, trans_type)

                        if (t < best_t) then
                           best_t = t
                           bi = c_bi; bj = c_bj; bk = c_bk; variant = c_var
                        end if
                     end do
                  end do
               end do
            else
               do ci = 1, size(n_bi_candidates_13)
                  c_bi = n_bi_candidates_13(ci)
                  if (c_bi > s1) cycle
                  do ck = 1, size(n_bk_candidates_13)
                     c_bk = n_bk_candidates_13(ck)
                     if (c_bk > s3) cycle
                     do cj = 1, size(n_bj_candidates_13)
                        c_bj = n_bj_candidates_13(cj)
                        if (c_bj > s2) cycle

                        ! L1 cache limit
                        if ( c_bi*c_bj*c_bk*ISFFT_ELEMSIZE() > 24576*3 ) cycle

                        t = Bench_Transpose_Once(ain, aout, s1,s2,s3, c_bi,c_bj,c_bk, c_var, trans_type)

                        if (t < best_t) then
                           best_t = t
                           bi = c_bi; bj = c_bj; bk = c_bk; variant = c_var
                        end if
                     end do
                  end do
               end do
            end if
         end do

         call Refine_Transpose_Around(ain, aout, s1,s2,s3, bi,bj,bk, variant, trans_type)

         ! Save results to corresponding variables
         if (trans_type == 12) then
            opt_lt12_bi = bi; opt_lt12_bj = bj; opt_lt12_bk = bk; opt_lt12_variant = variant
         else
            opt_lt13_bi = bi; opt_lt13_bj = bj; opt_lt13_bk = bk; opt_lt13_variant = variant
         end if
      end subroutine Tune_Transpose_BS
      !------------------------------------------------------------
      pure integer function ISFFT_ELEMSIZE() result(bytes)
         bytes = storage_size(0.0_ISFFT_REAL_KIND)/8
      end function ISFFT_ELEMSIZE
      !------------------------------------------------------------
      real(real64) function Bench_Transpose_Once(ain, aout, s1,s2,s3, bi,bj,bk, variant, trans_type) result(sec)
         real(kind=ISFFT_REAL_KIND), intent(in)  :: ain(:,:,:)
         real(kind=ISFFT_REAL_KIND), intent(out) :: aout(:,:,:)
         real(kind=ISFFT_REAL_KIND), allocatable :: temp(:,:,:)
         integer, intent(in) :: s1,s2,s3, bi,bj,bk, variant, trans_type
         integer(int64) :: c0,c1, rate
         integer :: r

         call Scrub_Transpose_Cache(ain, aout,  s1,s2,s3)

         if (trans_type == 12) then
            allocate(temp(s2,s1,s3))
         else
            allocate(temp(s3,s2,s1))
         end if

         call system_clock(count_rate=rate)
         call Kernel_Transpose_Test(ain, temp, s1,s2,s3, bi,bj,bk, variant, trans_type)
         call Kernel_Transpose_Test(temp, aout, s1,s2,s3, bi,bj,bk, variant, trans_type)
         call system_clock(c0)
         do r = 1, 2
            call Kernel_Transpose_Test(ain, temp, s1,s2,s3, bi,bj,bk, variant, trans_type)
            call Kernel_Transpose_Test(temp, aout, s1,s2,s3, bi,bj,bk, variant, trans_type)
         end do
         call system_clock(c1)

         sec = real(c1-c0,real64)/real(rate,real64) / 4.0_real64

         deallocate(temp)

      end function Bench_Transpose_Once
      !------------------------------------------------------------
      subroutine Kernel_Transpose_Test(ain, aout, dim1,dim2,dim3, bi,bj,bk, variant, trans_type)
         real(kind=ISFFT_REAL_KIND), intent(in)  :: ain(:,:,:)
         real(kind=ISFFT_REAL_KIND), intent(out) :: aout(:,:,:)
         integer, intent(in) :: dim1,dim2,dim3, bi,bj,bk, variant, trans_type
         integer :: ii, jj, kk, iend, jend, kend, i, j, k

         ! $omp parallel do schedule(static) collapse(3) private(ii,jj,kk,iend,jend,kend,i,j,k)
         do kk = 1, dim3, bk
            kend = min(kk+bk-1, dim3)
            do jj = 1, dim2, bj
               jend = min(jj+bj-1, dim2)
               do ii = 1, dim1, bi
                  iend = min(ii+bi-1, dim1)

                  if (trans_type == 12) then
                     if (variant == 0) then
                        ! 0：k→i→j, sequential write
                        do k = kk, kend
                           do i = ii, iend
                              !$omp simd
                              do j = jj, jend
                                 aout(j,i,k)   = ain(i,j,k)
                              end do
                           end do
                        end do
                     else
                        ! 1：k→j→i, sequential read
                        do k = kk, kend
                           do j = jj, jend
                              !$omp simd
                              do i = ii, iend
                                 aout(j,i,k)   = ain(i,j,k)
                              end do
                           end do
                        end do
                     end if
                  else  ! trans_type == 13
                     if (variant == 0) then
                        ! 0：i→j→k, sequential write
                        do i = ii, iend
                           do j = jj, jend
                              !$omp simd
                              do k = kk, kend
                                 aout(k,j,i)   = ain(i,j,k)
                              end do
                           end do
                        end do
                     else
                        ! 1：k→j→i, sequential read
                        do k = kk, kend
                           do j = jj, jend
                              !$omp simd
                              do i = ii, iend
                                 aout(k,j,i)   = ain(i,j,k)
                              end do
                           end do
                        end do
                     end if
                  end if

               end do
            end do
         end do
         ! $omp end parallel do
      end subroutine Kernel_Transpose_Test
      !------------------------------------------------------------
      subroutine Scrub_Transpose_Cache(a, b, s1,s2,s3)
         real(kind=ISFFT_REAL_KIND), intent(in)  :: a(:,:,:)
         real(kind=ISFFT_REAL_KIND), intent(out) :: b(:,:,:)
         integer, intent(in) :: s1,s2,s3
         integer :: i,j,k
         do k = 1, s3, max(1, s3/8)
            do j = 1, s2, max(1, s2/8)
               do i = 1, s1, max(1, s1/8)
                  b(1,1,1) = a(i,j,k)
               end do
            end do
         end do
      end subroutine Scrub_Transpose_Cache
      !------------------------------------------------------------
      subroutine Refine_Transpose_Around(ain, aout, s1,s2,s3, bi,bj,bk, variant, trans_type)
         real(kind=ISFFT_REAL_KIND), intent(in)  :: ain(:,:,:)
         real(kind=ISFFT_REAL_KIND), intent(out) :: aout(:,:,:)
         integer, intent(in) :: s1,s2,s3, variant, trans_type
         integer, intent(inout) :: bi,bj,bk
         integer :: c_bi, c_bj, c_bk
         real(real64) :: base_t, t

         base_t = Bench_Transpose_Once(ain, aout, s1,s2,s3, bi,bj,bk, variant, trans_type)

         ! Optimize different dimensions based on transpose type
         if (trans_type == 12) then
            ! For 12 transpose, optimize bi and bj
            do c_bi = max(8, bi-8), min(s1, bi+8), 8
               if (c_bi == bi) cycle
               t = Bench_Transpose_Once(ain, aout, s1,s2,s3, c_bi,bj,bk, variant, trans_type)
               if (t < base_t) then
                  bi = c_bi; base_t = t
               end if
            end do

            do c_bj = max(8, bj-8), min(s2, bj+8), 8
               if (c_bj == bj) cycle
               t = Bench_Transpose_Once(ain, aout, s1,s2,s3, bi,c_bj,bk, variant, trans_type)
               if (t < base_t) then
                  bj = c_bj; base_t = t
               end if
            end do
         else  ! trans_type == 13
            ! For 13 transpose, optimize bi and bk
            do c_bi = max(8, bi-8), min(s1, bi+8), 8
               if (c_bi == bi) cycle
               t = Bench_Transpose_Once(ain, aout, s1,s2,s3, c_bi,bj,bk, variant, trans_type)
               if (t < base_t) then
                  bi = c_bi; base_t = t
               end if
            end do

            do c_bk = max(8, bk-8), min(s3, bk+8), 8
               if (c_bk == bk) cycle
               t = Bench_Transpose_Once(ain, aout, s1,s2,s3, bi,bj,c_bk, variant, trans_type)
               if (t < base_t) then
                  bk = c_bk; base_t = t
               end if
            end do
         end if
      end subroutine Refine_Transpose_Around
   end subroutine Tune_LocalTranspose
   !----------------------------------------------------------------------
   subroutine Reset_LocalTranspose()
      implicit none

      enable_transpose_tuning = .true.
      If_lt12_tuned=.false.
      If_lt13_tuned=.false.


   end subroutine Reset_LocalTranspose
   !----------------------------------------------------------------------

end module ISFFT_LT_Tuning
