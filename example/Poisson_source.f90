!-------------------------------------------------------------------------------
subroutine Poisson_source(f)
   use CFD_solver_parameters
   implicit none

   real(kind=CFD_solver_real_kind), intent(out) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
   real(kind=CFD_solver_real_kind), parameter :: pi = 3.14159265358979323846_CFD_solver_real_kind
   real(kind=CFD_solver_real_kind) :: axis_mode(max(nx,ny,nz),3)
   real(kind=CFD_solver_real_kind) :: wave_number(3), domain_length(3), grid_spacing(3)
   real(kind=CFD_solver_real_kind) :: coordinate, laplacian_factor
   integer :: local_size(3), global_offset(3)
   integer :: axis, node, i, j, k

   local_size = [nx, ny, nz]
   global_offset = [i_offset(npx), j_offset(npy), k_offset(npz)]
   domain_length = [SLx, SLy, SLz]
   grid_spacing = [hx, hy, hz]
   axis_mode = 1.0_CFD_solver_real_kind
   wave_number = 0.0_CFD_solver_real_kind

   ! P uses cos(2*pi*x/L), D uses sin(pi*x/L), N uses cos(pi*x/L).
   ! Inactive directions contribute a factor of one and no derivative.
   do axis = 1, Dimension
      wave_number(axis) = pi/domain_length(axis)
      if (Boundary(axis) == 0) wave_number(axis) = 2.0_CFD_solver_real_kind*wave_number(axis)
      do node = 1, local_size(axis)
         coordinate = real(global_offset(axis)+node-2, CFD_solver_real_kind)*grid_spacing(axis)
         if (Boundary(axis) == 1) then
            axis_mode(node,axis) = sin(wave_number(axis)*coordinate)
         else
            axis_mode(node,axis) = cos(wave_number(axis)*coordinate)
         end if
      end do
   end do

   ! phi is the product of the directional modes, and f is its continuous Laplacian.
   laplacian_factor = -sum(wave_number**2)
   f = 0.0_CFD_solver_real_kind
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            f(i,j,k) = laplacian_factor*axis_mode(i,1)*axis_mode(j,2)*axis_mode(k,3)
         end do
      end do
   end do
end subroutine Poisson_source
!-------------------------------------------------------------------------------
