!------------------------------------------------------------------------------
! ISFFT_Redistribution
! Public redistribution and decomposition-activation entry points.
!------------------------------------------------------------------------------
module ISFFT_Redistribution
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   use ISFFT_Problem_State, only: internal_layout, original_layout
   use ISFFT_MPI_Part, only: ISFFT_part_change
   use ISFFT_Communication, only: Redistribute_Field, Require_Redistribution_Initialization
   implicit none
   private
   public :: ISFFT_Activate_Internal_Layout, ISFFT_Activate_Original_Layout, &
      ISFFT_Redistribution_Forward, ISFFT_Redistribution_Inverse

contains

   !------------------------------------------------------------------------------
   !! Public movement and layout-activation entry points
   !------------------------------------------------------------------------------
   subroutine ISFFT_Redistribution_Forward(field)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: field(:, :, :)

      call Redistribute_Field(field, .true.)
   end subroutine ! ISFFT_Redistribution_Forward

   subroutine ISFFT_Redistribution_Inverse(field)
      real(kind=ISFFT_REAL_KIND), allocatable, intent(inout) :: field(:, :, :)

      call Redistribute_Field(field, .false.)
   end subroutine ! ISFFT_Redistribution_Inverse

   !! Metadata-only activation. Caller storage must already match the target layout.
   subroutine ISFFT_Activate_Internal_Layout()

      call Require_Redistribution_Initialization()
      call ISFFT_part_change(internal_layout%process_grid(1), internal_layout%process_grid(2), &
         internal_layout%process_grid(3))
   end subroutine ! ISFFT_Activate_Internal_Layout

   !! Metadata-only activation. Caller storage must already match the original layout.
   subroutine ISFFT_Activate_Original_Layout()

      call Require_Redistribution_Initialization()
      call ISFFT_part_change(original_layout%process_grid(1), original_layout%process_grid(2), &
         original_layout%process_grid(3))
   end subroutine ! ISFFT_Activate_Original_Layout

end module ISFFT_Redistribution
