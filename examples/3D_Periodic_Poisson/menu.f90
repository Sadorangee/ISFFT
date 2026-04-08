!!--------------------------------------------------------------------------------
include'parameters.f90'
include'part_mpi.f90'
include'main_solver.f90'

Program menu
use flow_parameters
use TMS_constants
implicit none
integer :: ierr, status(MPI_status_size)
!---------------------------------------------------------------------------------
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
call mpi_comm_size(MPI_COMM_WORLD,np_size,ierr)
call read_parameters

!------ npx0,npy0,npz0 are parallel partitions in x, y and z directions----------
call part

call main_solver

call mpi_finalize(ierr)

end Program menu
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine read_parameters
   use flow_parameters
   implicit none
   integer :: i, j, k, nk, nr, ierr, ntmp(100)
   real(kind=TMS_REAL_KIND) :: rtmp(100)
   real(kind=TMS_REAL_KIND), parameter :: PI=3.14159265358979d0
   !-------------------------------------------
   nparameters=0
   rparameters=0.d0

   If(my_id .eq. 0) then

      open(30,file='Poisson_solver_UI.in')

      write(*,*) 'check1 ===='
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*) nx_global, ny_global, nz_global
      read(30,*)
      read(30,*) npx0,npy0,npz0,LAP
      read(30,*)
      read(30,*) SLx, SLy, SLz
      read(30,*)
      read(30,*) CFL, T
      read(30,*)
      read(30,*) If_3dfft_decomp
      read(30,*)
      read(30,*) If_scheme
      read(30,*)
      close(30)

   Endif


   ntmp(1)= nx_global
   ntmp(2)= ny_global
   ntmp(3)= nz_global
   ntmp(4)= LAP
   ntmp(5)= npx0
   ntmp(6)= npy0
   ntmp(7)= npz0
   ntmp(8)= If_3dfft_decomp
   ntmp(9)= If_scheme


   call MPI_bcast(ntmp(1),100,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)


   nx_global  = ntmp(1)
   ny_global  = ntmp(2)
   nz_global  = ntmp(3)
   LAP  = ntmp(4)
   npx0= ntmp(5)
   npy0= ntmp(6)
   npz0= ntmp(7)
   If_3dfft_decomp    = ntmp(8)
   If_scheme          = ntmp(9)


   rtmp(1)= SLx
   rtmp(2)= SLy
   rtmp(3)= SLz
   rtmp(4)= CFL
   rtmp(5)= T

   call MPI_bcast(rtmp(1),100,TMS_DATA_TYPE,0, MPI_COMM_WORLD,ierr)

   SLx= rtmp(1)
   SLy= rtmp(2)
   SLz= rtmp(3)
   CFL= rtmp(4)
   T= rtmp(5)

   if(my_id.eq.0) then

      print*, "Poisson Solver by Li Qing and Li Zecheng"
      print*, "Mesh (Nx,Ny,Nz): " , nx_global, ny_global, nz_global
      print*, "3D Partition: "    ,npx0,npy0,npz0, "  Total procs= ",npx0*npy0*npz0
      if (If_scheme .eq. 0) print*, "Scheme: Pseudo-Spectral"  
      if (If_scheme .eq. 1) print*, "Scheme: FDM (2nd CD)"  

      if (np_size .eq. 1) then
         print*, "Serial 3D FFT code"
      else if (np_size .gt. 1) then
         if (If_3dfft_decomp .eq. 0) print*, "Block for parallel 3D FFT code"
         if (If_3dfft_decomp .eq. 1 ) print*, "Slab for parallel 3D FFT code"
         if (If_3dfft_decomp .eq. 2 ) print*, "Pencil for parallel 3D FFT code"

      endif

   endif


end subroutine read_parameters




