!------------------------------------------------------------------------------
! ISFFT_Discrete_Operators
! Discrete-operator coefficients and modal factors.
!------------------------------------------------------------------------------
module ISFFT_Discrete_Operators
   use ISFFT_Precision, only: ISFFT_REAL_KIND
   implicit none
   private
   public :: ISFFT_BC_PERIODIC, ISFFT_BC_DIRICHLET, ISFFT_BC_NEUMANN
   public :: ISFFT_Scheme_Coefficients, ISFFT_Get_Scheme_Coefficients
   public :: ISFFT_Lambda_Periodic, ISFFT_Lambda_Nonperiodic, ISFFT_Lambda_For_BC

   integer, parameter :: ISFFT_BC_PERIODIC = 0
   integer, parameter :: ISFFT_BC_DIRICHLET = 1
   integer, parameter :: ISFFT_BC_NEUMANN = 2

   real(kind=ISFFT_REAL_KIND), parameter :: ISFFT_PI = 3.1415926535897932384626433832795_ISFFT_REAL_KIND

   type :: ISFFT_Scheme_Coefficients
      real(kind=ISFFT_REAL_KIND) :: alpha = 0.0_ISFFT_REAL_KIND
      real(kind=ISFFT_REAL_KIND) :: beta = 0.0_ISFFT_REAL_KIND
      real(kind=ISFFT_REAL_KIND) :: a = 0.0_ISFFT_REAL_KIND
      real(kind=ISFFT_REAL_KIND) :: b = 0.0_ISFFT_REAL_KIND
      real(kind=ISFFT_REAL_KIND) :: c = 0.0_ISFFT_REAL_KIND
   end type ISFFT_Scheme_Coefficients

contains

   !------------------------------------------------------------------------------
   !! Discretization coefficients
   !------------------------------------------------------------------------------
   subroutine ISFFT_Get_Scheme_Coefficients(discretization_scheme, coeff)
      integer, intent(in) :: discretization_scheme
      type(ISFFT_Scheme_Coefficients), intent(out) :: coeff

      select case (discretization_scheme)
      case (1)
         coeff%alpha = 0.0_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 1.0_ISFFT_REAL_KIND
         coeff%b = 0.0_ISFFT_REAL_KIND
         coeff%c = 0.0_ISFFT_REAL_KIND
      case (2)
         coeff%alpha = 0.0_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 4.0_ISFFT_REAL_KIND / 3.0_ISFFT_REAL_KIND
         coeff%b = -1.0_ISFFT_REAL_KIND / 3.0_ISFFT_REAL_KIND
         coeff%c = 0.0_ISFFT_REAL_KIND
      case (3)
         coeff%alpha = 0.1_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 6.0_ISFFT_REAL_KIND / 5.0_ISFFT_REAL_KIND
         coeff%b = 0.0_ISFFT_REAL_KIND
         coeff%c = 0.0_ISFFT_REAL_KIND
      case (4)
         coeff%alpha = 0.0_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 3.0_ISFFT_REAL_KIND / 2.0_ISFFT_REAL_KIND
         coeff%b = -3.0_ISFFT_REAL_KIND / 5.0_ISFFT_REAL_KIND
         coeff%c = 1.0_ISFFT_REAL_KIND / 10.0_ISFFT_REAL_KIND
      case (5)
         coeff%alpha = 2.0_ISFFT_REAL_KIND / 11.0_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 12.0_ISFFT_REAL_KIND / 11.0_ISFFT_REAL_KIND
         coeff%b = 3.0_ISFFT_REAL_KIND / 11.0_ISFFT_REAL_KIND
         coeff%c = 0.0_ISFFT_REAL_KIND
      case (6)
         !! Double-precision spectral-like coefficients. The shortened
         !! single-precision values impose an artificial error floor near 1e-8.
         coeff%alpha = 0.502092570454933_ISFFT_REAL_KIND
         coeff%beta = 0.055691658873965_ISFFT_REAL_KIND
         coeff%a = 0.215649631004306_ISFFT_REAL_KIND
         coeff%b = 1.723321721695239_ISFFT_REAL_KIND
         coeff%c = 0.176597105958252_ISFFT_REAL_KIND
      case default
         coeff%alpha = 0.0_ISFFT_REAL_KIND
         coeff%beta = 0.0_ISFFT_REAL_KIND
         coeff%a = 1.0_ISFFT_REAL_KIND
         coeff%b = 0.0_ISFFT_REAL_KIND
         coeff%c = 0.0_ISFFT_REAL_KIND
      end select
   end subroutine ISFFT_Get_Scheme_Coefficients

   !------------------------------------------------------------------------------
   !! Boundary-specific spectral symbols
   !------------------------------------------------------------------------------
   pure function ISFFT_Lambda_Periodic(mode_index, n_grid, length, discretization_scheme, coeff) result(lambda)
      integer, intent(in) :: mode_index, n_grid, discretization_scheme
      real(kind=ISFFT_REAL_KIND), intent(in) :: length
      type(ISFFT_Scheme_Coefficients), intent(in) :: coeff
      real(kind=ISFFT_REAL_KIND) :: lambda
      real(kind=ISFFT_REAL_KIND) :: theta, h

      if (discretization_scheme == 0) then
         lambda = (2.0_ISFFT_REAL_KIND * ISFFT_PI * real(mode_index, ISFFT_REAL_KIND) / length)**2
      else
         h = length / real(n_grid, ISFFT_REAL_KIND)
         theta = 2.0_ISFFT_REAL_KIND * ISFFT_PI * real(mode_index, ISFFT_REAL_KIND) / real(n_grid, ISFFT_REAL_KIND)
         lambda = (4.0_ISFFT_REAL_KIND * coeff%a * sin(0.5_ISFFT_REAL_KIND * theta)**2 + coeff%b * sin(theta)**2 + &
            (4.0_ISFFT_REAL_KIND * coeff%c / 9.0_ISFFT_REAL_KIND) * sin(1.5_ISFFT_REAL_KIND * theta)**2) / &
            (h**2 * (1.0_ISFFT_REAL_KIND + 2.0_ISFFT_REAL_KIND * coeff%alpha * cos(theta) + 2.0_ISFFT_REAL_KIND * &
            coeff%beta * cos(2.0_ISFFT_REAL_KIND * theta)))
      end if
   end function ISFFT_Lambda_Periodic

   pure function ISFFT_Lambda_Nonperiodic(mode_index, n_grid, length, discretization_scheme, coeff) result(lambda)
      integer, intent(in) :: mode_index, n_grid, discretization_scheme
      real(kind=ISFFT_REAL_KIND), intent(in) :: length
      type(ISFFT_Scheme_Coefficients), intent(in) :: coeff
      real(kind=ISFFT_REAL_KIND) :: lambda
      real(kind=ISFFT_REAL_KIND) :: theta, h

      if (n_grid <= 1) then
         lambda = 0.0_ISFFT_REAL_KIND
         return
      end if

      if (discretization_scheme == 0) then
         lambda = (ISFFT_PI * real(mode_index, ISFFT_REAL_KIND) / length)**2
      else
         h = length / real(n_grid - 1, ISFFT_REAL_KIND)
         theta = ISFFT_PI * real(mode_index, ISFFT_REAL_KIND) / real(n_grid - 1, ISFFT_REAL_KIND)
         lambda = (4.0_ISFFT_REAL_KIND * coeff%a * sin(0.5_ISFFT_REAL_KIND * theta)**2 + coeff%b * sin(theta)**2 + &
            (4.0_ISFFT_REAL_KIND * coeff%c / 9.0_ISFFT_REAL_KIND) * sin(1.5_ISFFT_REAL_KIND * theta)**2) / &
            (h**2 * (1.0_ISFFT_REAL_KIND + 2.0_ISFFT_REAL_KIND * coeff%alpha * cos(theta) + 2.0_ISFFT_REAL_KIND * &
            coeff%beta * cos(2.0_ISFFT_REAL_KIND * theta)))
      end if
   end function ISFFT_Lambda_Nonperiodic

   pure function ISFFT_Lambda_For_BC(boundary_type, mode_index, n_grid, length, discretization_scheme, coeff) &
      result(lambda)
      integer, intent(in) :: boundary_type, mode_index, n_grid, discretization_scheme
      real(kind=ISFFT_REAL_KIND), intent(in) :: length
      type(ISFFT_Scheme_Coefficients), intent(in) :: coeff
      real(kind=ISFFT_REAL_KIND) :: lambda

      select case (boundary_type)
      case (ISFFT_BC_PERIODIC)
         lambda = ISFFT_Lambda_Periodic(mode_index, n_grid, length, discretization_scheme, coeff)
      case (ISFFT_BC_DIRICHLET, ISFFT_BC_NEUMANN)
         lambda = ISFFT_Lambda_Nonperiodic(mode_index, n_grid, length, discretization_scheme, coeff)
      case default
         lambda = 0.0_ISFFT_REAL_KIND
      end select
   end function ISFFT_Lambda_For_BC

end module ISFFT_Discrete_Operators
