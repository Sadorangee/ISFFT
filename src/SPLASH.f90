module SPLASH
   use SPLASH_Planner, only: SPLASH_Initialize, SPLASH_Finalize
   use SPLASH_Integrated_Solver, only: SPLASH_Periodic_Poisson_solver
   implicit none
   
   private
   public :: SPLASH_Initialize, SPLASH_Periodic_Poisson_solver, SPLASH_Finalize

end module SPLASH

