!------------------------------------------------------------------------------
! ISFFT_Integrated_Solver
! Integrated equation-solver workflows.
!------------------------------------------------------------------------------
module ISFFT_Integrated_Solver
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   use ISFFT_Plans, only: poisson_plan, redistribution_forward_plan, redistribution_inverse_plan
   use ISFFT_Problem_State, only: boundary_types, dimensions, halo_width, internal_layout, &
      npx0, npy0, npz0, nx, &
      ny, nz, original_layout
   use ISFFT_Workspace, only: Load_Field, poisson_field, poisson_symbol
   use ISFFT_Halos_Exchange, only: ISFFT_Update_Boundary_x, ISFFT_Update_Boundary_xy, ISFFT_Update_Boundary_xyz
   use ISFFT_Modal_Operations, only: Apply_Poisson_Operator
   use ISFFT_Communication, only: Move_State
   use ISFFT_Transform_Execution, only: Can_Direct_Input, Can_Direct_Output, Execute_Forward, Execute_Inverse
   implicit none
   private
   public :: ISFFT_Poisson_solver

contains

   !------------------------------------------------------------------------------
   !! Poisson solve in one, two or three dimensions
   !------------------------------------------------------------------------------
   subroutine ISFFT_Poisson_solver(phi, f)
      real(kind=ISFFT_REAL_KIND), intent(inout) &
         :: phi(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)
      real(kind=ISFFT_REAL_KIND), intent(in) &
         :: f(1 - halo_width:nx + halo_width, 1 - halo_width:ny + halo_width, 1 - halo_width:nz + halo_width)
      logical :: needs_redistribution, direct_output, direct_input

      if (.not. poisson_plan%ready) error stop 'ISFFT: call ISFFT_Initialize before Poisson solve'
      if (dimensions < 1 .or. dimensions > 3) error stop 'invalid Poisson dimensions'
      needs_redistribution = any([npx0, npy0, npz0] /= internal_layout%process_grid)
      if (needs_redistribution .and. any([npx0, npy0, npz0] /= original_layout%process_grid)) &
         error stop 'Poisson input layout is neither caller nor target layout'
      associate (field => poisson_field)

         !! Forward transform from caller data or a reusable field buffer.
         direct_input = halo_width == 0 .and. .not. needs_redistribution .and. Can_Direct_Input(poisson_plan)
         if (direct_input) then
            call Execute_Forward(poisson_plan, field, f)
         else
            call Load_Field(field, f(1:nx, 1:ny, 1:nz))
            if (needs_redistribution) call Move_State(redistribution_forward_plan, field)
            call Execute_Forward(poisson_plan, field)
         end if

         call Apply_Poisson_Operator(poisson_symbol, field)

         !! Inverse transform and restore caller ownership.
         direct_output = halo_width == 0 .and. .not. needs_redistribution .and. Can_Direct_Output(poisson_plan)
         if (direct_output) then
            call Execute_Inverse(poisson_plan, field, phi)
         else
            call Execute_Inverse(poisson_plan, field)
            if (needs_redistribution) call Move_State(redistribution_inverse_plan, field)
            phi = 0
            phi(1:nx, 1:ny, 1:nz) = field%r
         end if

         select case (dimensions)
         case (1)
            call ISFFT_Update_Boundary_x(phi, boundary_types(1))
         case (2)
            call ISFFT_Update_Boundary_xy(phi)
         case (3)
            call ISFFT_Update_Boundary_xyz(phi)
         end select
      end associate
   end subroutine ! ISFFT_Poisson_solver

end module ISFFT_Integrated_Solver
