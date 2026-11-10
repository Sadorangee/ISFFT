!------------------------------------------------------------------------------
! ISFFT_Modal_Operations
! Pointwise modal operations using prepared coefficients.
!------------------------------------------------------------------------------
module ISFFT_Modal_Operations
   use ISFFT_Workspace, only: Poisson_Operator, Spectral_Field
   implicit none
   private
   public :: Apply_Poisson_Operator

contains

   subroutine Apply_Poisson_Operator(op, field)
      type(Poisson_Operator), intent(in) :: op
      type(Spectral_Field), intent(inout) :: field

      if (allocated(field%z)) then
         field%z = field%z * op%multiplier
      else
         field%r = field%r * op%multiplier
      end if
   end subroutine ! Apply_Poisson_Operator

end module ISFFT_Modal_Operations
