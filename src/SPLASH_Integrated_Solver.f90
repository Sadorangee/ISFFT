!!----------------------------------------------------------------------
module SPLASH_Integrated_Solver
   use SPLASH_Precision
   use SPLASH_Parameters
   use SPLASH_MPI_Constants
   use SPLASH_Buffer
   use SPLASH_Repartitioning
   use SPLASH_3D_Periodic_Poisson
   use SPLASH_Halos_Exchange
   implicit none

   private
   public :: SPLASH_Periodic_Poisson_solver

contains

!!----------------------------------------------------------------------
   subroutine SPLASH_Periodic_Poisson_solver(phi, f)
      implicit none

      real(kind=SPLASH_REAL_KIND),dimension(:,:,:), intent(inout) :: phi, f
      !!!---------------------------------------------------------------
      allocate(f_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
      f_internal = f

      call SPLASH_Repartitioning_Forward(f_internal)

      allocate(phi_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))

      if ( np_size .eq. 1 ) call SPLASH_3D_Periodic_Poisson_serial(phi_internal,f_internal)

      if ( np_size .gt. 1 ) then

         if (If_3dfft_decomp .eq. 0) call SPLASH_3D_Periodic_Poisson_block(phi_internal,f_internal)
         if (If_3dfft_decomp .eq. 1) call SPLASH_3D_Periodic_Poisson_slab(phi_internal,f_internal)
         if (If_3dfft_decomp .eq. 2) call SPLASH_3D_Periodic_Poisson_pencil(phi_internal,f_internal)

      end if

      call SPLASH_Repartitioning_Inverse(phi_internal)

      call SPLASH_Update_Periodic_Boundary_xyz(phi_internal)

      phi = phi_internal
      deallocate(f_internal, phi_internal)

   End subroutine SPLASH_Periodic_Poisson_solver
!!-----------------------------------------------------------------------

end module SPLASH_Integrated_Solver
!!-----------------------------------------------------------------------
