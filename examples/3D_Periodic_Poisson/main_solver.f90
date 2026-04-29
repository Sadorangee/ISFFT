!!!----------------------------------------
include'Poisson_3D_source.f90'
include'IO_Poisson_3D.f90'
!!!----------------------------------------
subroutine main_solver
   Use flow_parameters
   use TMS_constants
   use ISFFT
   implicit none

   real(kind=TMS_REAL_KIND), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP):: phi,f
   integer :: Global_Grid_Size(3)
   integer :: Computational_Domain_Size(3)
   integer :: Process_Grid_Size(3)
   integer :: ierr

   phi = 0.d0
   f = 0.d0

   hx=SLx/nx_global
   hy=SLy/ny_global
   hz=SLz/nz_global

   dt = CFL*min(hx,hy,hz)**2.0
   Nt = max(1, int(T/dt))

   Global_Grid_Size = [nx_global, ny_global, nz_global]
   Computational_Domain_Size = [SLx, SLy, SLz]
   Process_Grid_Size = [npx0, npy0, npz0]
   
   Nt = 1 !! Manually set the number of loops. 

   call ISFFT_Initialize(Global_Grid_Size, Process_Grid_Size, Computational_Domain_Size, &
      LAP, If_3dfft_decomp, If_scheme, my_id)

   !!!---------------------------------------------------------------
   !!!-----------------------Main solver-----------------------------
   !!!---------------------------------------------------------------
   do Istep = 1, Nt

      ! In real CFD simulations, the source term f is updated at each time step.
      call Poisson_3D_source(f)
      call ISFFT_Periodic_Poisson_solver(phi, f)

   end do
   !!!---------------------------------------------------------------
   !!!---------------------------------------------------------------
   !!!---------------------------------------------------------------
   call ISFFT_Finalize()

   call IO_Poisson_3D(phi)
   if (my_id .eq. 0) write(*,*) '3D Periodic Poisson Solve is done!!'


End subroutine main_solver
!!!----------------------------------------
