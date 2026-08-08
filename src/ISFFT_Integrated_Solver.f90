!!----------------------------------------------------------------------
module ISFFT_Integrated_Solver
   use ISFFT_Precision
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_Repartitioning
   use ISFFT_1D_Mixed_Poisson
   use ISFFT_2D_Mixed_Poisson
   use ISFFT_3D_Mixed_Poisson
   use ISFFT_Halos_Exchange
   implicit none

   private
   public :: ISFFT_Poisson_solver, ISFFT_1D_Poisson_solver, ISFFT_2D_Poisson_solver
   public :: ISFFT_3D_Poisson_solver

contains

!!----------------------------------------------------------------------
   subroutine ISFFT_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      select case (Dimension)
      case (1)
         call ISFFT_1D_Unified_Poisson_solver(phi, f)
      case (2)
         call ISFFT_2D_Unified_Poisson_solver(phi, f)
      case (3)
         call ISFFT_3D_Unified_Poisson_solver(phi, f)
      end select
   end subroutine ISFFT_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_1D_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      call ISFFT_1D_Unified_Poisson_solver(phi, f)
   end subroutine ISFFT_1D_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_1D_Unified_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      allocate(f_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
      f_internal = 0.0_ISFFT_REAL_KIND
      f_internal(1:nx,1,1) = f(1:nx,1,1)

      allocate(phi_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
         source=0.0_ISFFT_REAL_KIND)

      call ISFFT_1D_Mixed_Poisson_serial(phi_internal, f_internal)
      call ISFFT_Update_Boundary_x(phi_internal, BC(1))

      phi = phi_internal
      deallocate(f_internal, phi_internal)
   end subroutine ISFFT_1D_Unified_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_2D_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      call ISFFT_2D_Unified_Poisson_solver(phi, f)
   end subroutine ISFFT_2D_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_2D_Unified_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      allocate(f_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), source=0.0_ISFFT_REAL_KIND)
      f_internal(1:nx,1:ny,1) = f(1:nx,1:ny,1)
      call ISFFT_Repartitioning_Forward(f_internal)

      allocate(phi_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), source=0.0_ISFFT_REAL_KIND)
      if (np_size == 1) then
         call ISFFT_2D_Mixed_Poisson_serial(phi_internal, f_internal)
      else if (If_decomp == 0) then
         call ISFFT_2D_Mixed_Poisson_block(phi_internal, f_internal)
      else
         call ISFFT_2D_Mixed_Poisson_slab(phi_internal, f_internal)
      end if

      call ISFFT_Repartitioning_Inverse(phi_internal)
      call ISFFT_Update_Boundary_xy(phi_internal)
      phi = phi_internal
      deallocate(f_internal, phi_internal)
   end subroutine ISFFT_2D_Unified_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_3D_Poisson_solver(phi, f)
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      call ISFFT_3D_Unified_Poisson_solver(phi, f)
   end subroutine ISFFT_3D_Poisson_solver

!!----------------------------------------------------------------------
   subroutine ISFFT_3D_Unified_Poisson_solver(phi, f)
      implicit none
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

      allocate(f_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
      f_internal = f

      call ISFFT_Repartitioning_Forward(f_internal)

      allocate(phi_internal(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
         source=0.0_ISFFT_REAL_KIND)

      if (np_size .eq. 1) then
         call ISFFT_3D_Mixed_Poisson_serial(phi_internal, f_internal)
      else if (If_decomp .eq. 0) then
         call ISFFT_3D_Mixed_Poisson_block(phi_internal, f_internal)
      else
         call ISFFT_3D_Mixed_Poisson_transpose(phi_internal, f_internal)
      end if

      call ISFFT_Repartitioning_Inverse(phi_internal)
      call ISFFT_Update_Boundary_xyz(phi_internal)

      phi = phi_internal
      deallocate(f_internal, phi_internal)
   end subroutine ISFFT_3D_Unified_Poisson_solver
!!----------------------------------------------------------------------

end module ISFFT_Integrated_Solver
!!----------------------------------------------------------------------
