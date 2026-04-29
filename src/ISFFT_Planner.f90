module ISFFT_Planner
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_MPI_Part, only: ISFFT_part
   use ISFFT_FFT_Pre, only: FFT_Initialize, FFT_Finalize
   use ISFFT_LT_Tuning, only: Tune_LocalTranspose, Reset_LocalTranspose
   implicit none

contains

   !----------------------------------------------------------------------
   subroutine ISFFT_Initialize(Global_Grid_Size, Process_Grid_Size, Computational_Domain_Size, &
      Num_Lap, If_Decomp_Type, If_Algebra_scheme, My_Rank)
      implicit none
      integer, intent(in) :: Global_Grid_Size(3)
      integer, intent(in) :: Process_Grid_Size(3)
      integer, intent(in) :: Computational_Domain_Size(3)
      integer, intent(in) :: Num_Lap
      integer, intent(in) :: If_Decomp_Type
      integer, intent(in) :: If_Algebra_scheme
      integer, intent(in) :: My_Rank
      integer :: ierr

      !! ------------------Read parameters--------------------
      nx_global = Global_Grid_Size(1)
      ny_global = Global_Grid_Size(2)
      nz_global = Global_Grid_Size(3)

      npx0 = Process_Grid_Size(1)
      npy0 = Process_Grid_Size(2)
      npz0 = Process_Grid_Size(3)

      SLx = Computational_Domain_Size(1)
      SLy = Computational_Domain_Size(2)
      SLz = Computational_Domain_Size(3)

      LAP = Num_Lap

      If_3dfft_decomp = If_Decomp_Type
      If_scheme = If_Algebra_scheme

      my_id = My_Rank
      np_size = npx0 * npy0 * npz0

      halfspectrum_x = nx_global / 2 + 1
      halfspectrum_y = ny_global / 2 + 1
      halfspectrum_z = nz_global / 2 + 1

      call ISFFT_part()

      if (np_size .eq. 1 .or. If_3dfft_decomp .eq. 0) enable_transpose_tuning = .false.

      !! -----------------Initialize Buffer-------------------
      call Allocate_Buffer_3D_FFT()
      call Initialize_Buffer_3D_FFT()

      !! ------------------Initialize FFT---------------------
      call FFT_Initialize(nx_global, ny_global, nz_global)

      ! ---------------Tune Local Transpose------------------
      if (enable_transpose_tuning) then
         call Tune_LocalTranspose(12)
         call Tune_LocalTranspose(13)
      end if

      !! -----------------Print-------------------
      if (my_id .eq. 0) then
         write(*,*) 'ISFFT_Initialize is done!'
      end if

   end subroutine ISFFT_Initialize
   !----------------------------------------------------------------------
   subroutine ISFFT_Finalize()
      implicit none

      call Deallocate_Buffer_3D_FFT()

      call Reset_LocalTranspose()

      call FFT_Finalize()

      if (my_id .eq. 0) then
         write(*,*) 'ISFFT_Finalize is done!'
      end if

   end subroutine ISFFT_Finalize
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine Allocate_Buffer_3D_FFT()
      implicit none
      integer :: npx0_target, npy0_target, npz0_target
      integer :: nx_target, ny_target, nz_target
      integer :: ierr

      select case(If_3dfft_decomp)
       case (0)  ! Block Decomposition
         npx0_target = npx0
         npy0_target = npy0
         npz0_target = npz0

         nx_target = nx_global/npx0_target
         ny_target = ny_global/npy0_target
         nz_target = nz_global/npz0_target

       case (1)  ! Slab Decomposition
         npx0_target = 1
         npy0_target = 1
         npz0_target = npx0*npy0*npz0

         nx_target = nx_global/npx0_target
         ny_target = ny_global/npy0_target
         nz_target = nz_global/npz0_target

         if (np_size .gt. min(nx_global, ny_global, nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Slab Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

       case (2)  ! Pencil Decomposition
         npx0_target = 1
         npy0_target = 2**((int(log(real(npx0*npy0*npz0))/log(2.0)) + 1) / 2)
         npz0_target = npx0*npy0*npz0/npy0_target

         nx_target = nx_global/npx0_target
         ny_target = ny_global/npy0_target
         nz_target = nz_global/npz0_target

         if (np_size .gt. min(nx_global*ny_global, ny_global*nz_global, nx_global*nz_global)) then
            if(my_id .eq. 0) write(*,*) "Error: The Pencil Decomposition is not suitable."
            call mpi_finalize(ierr)
            stop
         endif

      end select

      if (Buffer_3DFFT_Allocated) return

      !! General FFT Buffer
      allocate(f_re(1:nx_target, 1:ny_target, 1:nz_target))
      allocate(f_im(1:nx_target, 1:ny_target, 1:nz_target))

      if (If_3dfft_decomp .ne. 0) then
         allocate(sendbuf_3DFFT(1:nx_target*ny_target*nz_target))
         allocate(sendbuf_im_3DFFT(1:nx_target*ny_target*nz_target))
         allocate(recvbuf_3DFFT(1:nx_target*ny_target*nz_target))
         allocate(recvbuf_im_3DFFT(1:nx_target*ny_target*nz_target))

         allocate(sendbuf_3DFFT_complex(1:2,1:nx_target*ny_target*nz_target))
         allocate(recvbuf_3DFFT_complex(1:2,1:nx_target*ny_target*nz_target))
      endif

      !! Block FFT Buffer
      if (If_3dfft_decomp .eq. 0) then

         allocate(x_sendbuf(1:nx_global))
         allocate(x_im_sendbuf(1:nx_global))
         allocate(x_recvbuf(1:nx_global))
         allocate(x_im_recvbuf(1:nx_global))

         allocate(y_sendbuf(1:ny_global))
         allocate(y_im_sendbuf(1:ny_global))
         allocate(y_recvbuf(1:ny_global))
         allocate(y_im_recvbuf(1:ny_global))

         allocate(z_sendbuf(1:nz_global))
         allocate(z_im_sendbuf(1:nz_global))
         allocate(z_recvbuf(1:nz_global))
         allocate(z_im_recvbuf(1:nz_global))

         !! Slab FFT Buffer
      else if (If_3dfft_decomp .eq. 1) then


         ! allocate(f_hat_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         ! allocate(f_hat_im_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(f_hat_y(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(f_hat_im_y(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(f_hat(1:nz_global/npy0_target, 1:nx_global/npx0_target, 1:ny_global/npz0_target))
         allocate(f_hat_im(1:nz_global/npy0_target, 1:nx_global/npx0_target, 1:ny_global/npz0_target))
         ! allocate(phi_hat_y(1:nz_global/npz0_target, 1:nx_global/npx0_target, 1:nz_global/npy0_target))
         ! allocate(phi_hat_im_y(1:nz_global/npz0_target, 1:nx_global/npx0_target, 1:nz_global/npy0_target))
         allocate(phi_hat_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(phi_hat_im_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(phi_re(1:nx_global/npx0_target, 1:ny_global/npy0_target, 1:nz_global/npz0_target))
         allocate(phi_im(1:nx_global/npx0_target, 1:ny_global/npy0_target, 1:nz_global/npz0_target))

         allocate(phi_re_transposed_a(1:nz_global/npz0_target, 1:nx_global/npx0_target, 1:ny_global/npy0_target))
         allocate(phi_im_transposed_a(1:nz_global/npz0_target, 1:nx_global/npx0_target, 1:ny_global/npy0_target))


      else if (If_3dfft_decomp .eq. 2) then


         ! allocate(f_hat_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         ! allocate(f_hat_im_x(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(f_hat_y(1:ny_global/npx0_target, 1:nx_global/npy0_target, 1:nz_global/npz0_target))
         allocate(f_hat_im_y(1:ny_global/npx0_target, 1:nx_global/npy0_target, 1:nz_global/npz0_target))
         allocate(f_hat(1:nz_global/npx0_target, 1:nx_global/npy0_target, 1:ny_global/npz0_target))
         allocate(f_hat_im(1:nz_global/npx0_target, 1:nx_global/npy0_target, 1:ny_global/npz0_target))

         ! allocate(phi_hat_y(1:nz_global/npz0_target, 1:nx_global/npy0_target, 1:ny_global/npx0_target))
         ! allocate(phi_hat_im_y(1:nz_global/npz0, 1:nx_global/npy0, 1:ny_global/npx0))
         allocate(phi_hat_x(1:ny_global/npx0_target, 1:nx_global/npy0_target, 1:nz_global/npz0_target))
         allocate(phi_hat_im_x(1:ny_global/npx0_target, 1:nx_global/npy0_target, 1:nz_global/npz0_target))
         allocate(phi_re(1:nx_global/npx0_target, 1:ny_global/npy0_target, 1:nz_global/npz0_target))
         allocate(phi_im(1:nx_global/npx0_target, 1:ny_global/npy0_target, 1:nz_global/npz0_target))

         allocate(phi_re_transposed_a(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(phi_im_transposed_a(1:ny_global/npy0_target, 1:nx_global/npx0_target, 1:nz_global/npz0_target))
         allocate(phi_re_transposed_b(1:nz_global/npz0_target, 1:nx_global/npy0_target, 1:ny_global/npx0_target))
         allocate(phi_im_transposed_b(1:nz_global/npz0_target, 1:nx_global/npy0_target, 1:ny_global/npx0_target))



      end if

      Buffer_3DFFT_Allocated = .true.

   end subroutine Allocate_Buffer_3D_FFT
   !----------------------------------------------------------------------
   subroutine Initialize_Buffer_3D_FFT()
      implicit none

      if (.not. Buffer_3DFFT_Allocated) return

      ! if (allocated(f_re)) f_re = 0.d0
      ! if (allocated(f_im)) f_im = 0.d0

      if (If_3dfft_decomp .eq. 0) then

         ! if (allocated(x_sendbuf)) x_sendbuf = 0.d0
         ! if (allocated(x_im_sendbuf)) x_im_sendbuf = 0.d0
         ! if (allocated(x_recvbuf)) x_recvbuf = 0.d0
         ! if (allocated(x_im_recvbuf)) x_im_recvbuf = 0.d0
         ! if (allocated(y_sendbuf)) y_sendbuf = 0.d0
         ! if (allocated(y_im_sendbuf)) y_im_sendbuf = 0.d0
         ! if (allocated(y_recvbuf)) y_recvbuf = 0.d0
         ! if (allocated(y_im_recvbuf)) y_im_recvbuf = 0.d0
         ! if (allocated(z_sendbuf)) z_sendbuf = 0.d0
         ! if (allocated(z_im_sendbuf)) z_im_sendbuf = 0.d0
         ! if (allocated(z_recvbuf)) z_recvbuf = 0.d0
         ! if (allocated(z_im_recvbuf)) z_im_recvbuf = 0.d0

      else if (If_3dfft_decomp .eq. 1) then

         ! if (allocated(sendbuf_3DFFT)) sendbuf_3DFFT = 0.d0
         ! if (allocated(sendbuf_im_3DFFT)) sendbuf_im_3DFFT = 0.d0
         ! if (allocated(recvbuf_3DFFT)) recvbuf_3DFFT = 0.d0
         ! if (allocated(recvbuf_im_3DFFT)) recvbuf_im_3DFFT = 0.d0


         ! if (allocated(f_hat_y)) f_hat_y = 0.d0
         ! if (allocated(f_hat_im_y)) f_hat_im_y = 0.d0
         ! if (allocated(f_hat)) f_hat = 0.d0
         ! if (allocated(f_hat_im)) f_hat_im = 0.d0
         ! if (allocated(phi_hat_z)) phi_hat_z = 0.d0
         ! if (allocated(phi_hat_im_z)) phi_hat_im_z = 0.d0
         ! if (allocated(phi_hat)) phi_hat = 0.d0
         ! if (allocated(phi_hat_im)) phi_hat_im = 0.d0

         ! if (allocated(phi_re_transposed_a)) phi_re_transposed_a = 0.d0
         ! if (allocated(phi_im_transposed_a)) phi_im_transposed_a = 0.d0



      else if (If_3dfft_decomp .eq. 2) then

         ! if (allocated(sendbuf_3DFFT)) sendbuf_3DFFT = 0.d0
         ! if (allocated(sendbuf_im_3DFFT)) sendbuf_im_3DFFT = 0.d0
         ! if (allocated(recvbuf_3DFFT)) recvbuf_3DFFT = 0.d0
         ! if (allocated(recvbuf_im_3DFFT)) recvbuf_im_3DFFT = 0.d0


         ! if (allocated(f_hat_x)) f_hat_x = 0.d0
         ! if (allocated(f_hat_im_x)) f_hat_im_x = 0.d0
         ! if (allocated(f_hat_y)) f_hat_y = 0.d0
         ! if (allocated(f_hat_im_y)) f_hat_im_y = 0.d0
         ! if (allocated(f_hat)) f_hat = 0.d0
         ! if (allocated(f_hat_im)) f_hat_im = 0.d0
         ! if (allocated(phi_hat_z)) phi_hat_z = 0.d0
         ! if (allocated(phi_hat_im_z)) phi_hat_im_z = 0.d0
         ! if (allocated(phi_hat_y)) phi_hat_y = 0.d0
         ! if (allocated(phi_hat_im_y)) phi_hat_im_y = 0.d0
         ! if (allocated(phi_hat)) phi_hat = 0.d0
         ! if (allocated(phi_hat_im)) phi_hat_im = 0.d0

         ! if (allocated(phi_re_transposed_a)) phi_re_transposed_a = 0.d0
         ! if (allocated(phi_im_transposed_a)) phi_im_transposed_a = 0.d0



      end if

   end subroutine Initialize_Buffer_3D_FFT
   !----------------------------------------------------------------------
   subroutine Deallocate_Buffer_3D_FFT()
      implicit none

      if (.not. Buffer_3DFFT_Allocated) return

      if (allocated(f_re)) deallocate(f_re)
      if (allocated(f_im)) deallocate(f_im)

      if (If_3dfft_decomp .ne. 0) then
         if (allocated(sendbuf_3DFFT)) deallocate(sendbuf_3DFFT)
         if (allocated(sendbuf_im_3DFFT)) deallocate(sendbuf_im_3DFFT)
         if (allocated(recvbuf_3DFFT)) deallocate(recvbuf_3DFFT)
         if (allocated(recvbuf_im_3DFFT)) deallocate(recvbuf_im_3DFFT)

         if (allocated(sendbuf_3DFFT_complex)) deallocate(sendbuf_3DFFT_complex)
         if (allocated(recvbuf_3DFFT_complex)) deallocate(recvbuf_3DFFT_complex)
      endif

      !! Block FFT Buffer
      if (If_3dfft_decomp .eq. 0) then

         if (allocated(x_sendbuf)) deallocate(x_sendbuf)
         if (allocated(x_im_sendbuf)) deallocate(x_im_sendbuf)
         if (allocated(x_recvbuf)) deallocate(x_recvbuf)
         if (allocated(x_im_recvbuf)) deallocate(x_im_recvbuf)

         if (allocated(y_sendbuf)) deallocate(y_sendbuf)
         if (allocated(y_im_sendbuf)) deallocate(y_im_sendbuf)
         if (allocated(y_recvbuf)) deallocate(y_recvbuf)
         if (allocated(y_im_recvbuf)) deallocate(y_im_recvbuf)

         if (allocated(z_sendbuf)) deallocate(z_sendbuf)
         if (allocated(z_im_sendbuf)) deallocate(z_im_sendbuf)
         if (allocated(z_recvbuf)) deallocate(z_recvbuf)
         if (allocated(z_im_recvbuf)) deallocate(z_im_recvbuf)

         !! Slab FFT Buffer
      else if (If_3dfft_decomp .eq. 1) then


         ! if (allocated(f_hat_x)) deallocate(f_hat_x)
         ! if (allocated(f_hat_im_x)) deallocate(f_hat_im_x)
         if (allocated(f_hat_y)) deallocate(f_hat_y)
         if (allocated(f_hat_im_y)) deallocate(f_hat_im_y)
         if (allocated(f_hat)) deallocate(f_hat)
         if (allocated(f_hat_im)) deallocate(f_hat_im)
         ! if (allocated(phi_hat_y)) deallocate(phi_hat_y)
         ! if (allocated(phi_hat_im_y)) deallocate(phi_hat_im_y)
         if (allocated(phi_hat_x)) deallocate(phi_hat_x)
         if (allocated(phi_hat_im_x)) deallocate(phi_hat_im_x)
         if (allocated(phi_re)) deallocate(phi_re)
         if (allocated(phi_im)) deallocate(phi_im)

         if (allocated(phi_re_transposed_a)) deallocate(phi_re_transposed_a)
         if (allocated(phi_im_transposed_a)) deallocate(phi_im_transposed_a)



      else if (If_3dfft_decomp .eq. 2) then


         if (allocated(f_hat_x)) deallocate(f_hat_x)
         if (allocated(f_hat_im_x)) deallocate(f_hat_im_x)
         if (allocated(f_hat_y)) deallocate(f_hat_y)
         if (allocated(f_hat_im_y)) deallocate(f_hat_im_y)
         if (allocated(f_hat)) deallocate(f_hat)
         if (allocated(f_hat_im)) deallocate(f_hat_im)

         if (allocated(phi_hat_y)) deallocate(phi_hat_y)
         if (allocated(phi_hat_im_y)) deallocate(phi_hat_im_y)
         if (allocated(phi_hat_x)) deallocate(phi_hat_x)
         if (allocated(phi_hat_im_x)) deallocate(phi_hat_im_x)
         if (allocated(phi_re)) deallocate(phi_re)
         if (allocated(phi_im)) deallocate(phi_im)

         if (allocated(phi_re_transposed_a)) deallocate(phi_re_transposed_a)
         if (allocated(phi_im_transposed_a)) deallocate(phi_im_transposed_a)
         if (allocated(phi_re_transposed_b)) deallocate(phi_re_transposed_b)
         if (allocated(phi_im_transposed_b)) deallocate(phi_im_transposed_b)



      end if

      Buffer_3DFFT_Allocated = .false.

   end subroutine Deallocate_Buffer_3D_FFT
   !----------------------------------------------------------------------



end module ISFFT_Planner
