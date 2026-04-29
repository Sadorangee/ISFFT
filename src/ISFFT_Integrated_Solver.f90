!!----------------------------------------------------------------------
module ISFFT_Integrated_Solver
   use ISFFT_Precision
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_Repartitioning
   use ISFFT_3D_Periodic_Poisson
   use ISFFT_Halos_Exchange
   implicit none

   private
   public :: ISFFT_Periodic_Poisson_solver

contains

!!----------------------------------------------------------------------
   subroutine ISFFT_Periodic_Poisson_solver(phi, f)
      implicit none

      ! real(kind=ISFFT_REAL_KIND),dimension(:,:,:), intent(inout) :: phi, f
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      !!!---------------------------------------------------------------
      allocate(f_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
      f_internal = f

      call ISFFT_Repartitioning_Forward(f_internal)

      allocate(phi_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))

      if ( np_size .eq. 1 ) call ISFFT_3D_Periodic_Poisson_serial(phi_internal,f_internal)

      if ( np_size .gt. 1 ) then

         if (If_3dfft_decomp .eq. 0) call ISFFT_3D_Periodic_Poisson_block(phi_internal,f_internal)
         if (If_3dfft_decomp .eq. 1) call ISFFT_3D_Periodic_Poisson_slab(phi_internal,f_internal)
         if (If_3dfft_decomp .eq. 2) call ISFFT_3D_Periodic_Poisson_pencil(phi_internal,f_internal)

      end if

      call ISFFT_Repartitioning_Inverse(phi_internal)

      call ISFFT_Update_Periodic_Boundary_xyz(phi_internal)

      phi = phi_internal
      deallocate(f_internal, phi_internal)

   End subroutine ISFFT_Periodic_Poisson_solver
!!-----------------------------------------------------------------------

end module ISFFT_Integrated_Solver
!!-----------------------------------------------------------------------
