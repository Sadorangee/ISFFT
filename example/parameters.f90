!!---------------------------------------------------------------------------------------------------------
module CFD_solver_precision
   use mpi
   implicit none
!------For Doubleprecision  (real*8)------------------------------------------------------------------
   integer,parameter::CFD_solver_real_kind=8,  CFD_solver_real_mpi_type=MPI_DOUBLE_PRECISION   ! double precison computing
   integer,parameter::CFD_solver_integer_mpi_type=MPI_INTEGER
! ------For Single precision (real*4)-----------------------
   ! integer,parameter::CFD_solver_real_kind=4,  CFD_solver_real_mpi_type=MPI_REAL             !  single precision computing
!===========  Parameters for MPI ==========================================================----------------
end module CFD_solver_precision

!!------parameters used in the CFD toy solver---------------------------------------
module CFD_solver_constants
!!---------------------------------------------------------------------------------
   Use CFD_solver_precision
   implicit none
!!-------------------------------------------------------------------------
   real(kind=CFD_solver_real_kind),parameter:: TOL_matrix=1e-6;
!--------MPI-------------------
   integer,save:: my_id,npx,npy,npz,npx0,npy0,npz0, ID_XP1,ID_XM1,ID_YP1,ID_YM1,ID_ZP1,ID_ZM1, &
       CFD_solver_barrier_level,np_size
   integer,save::  i_offset(0:2048),j_offset(0:2048),k_offset(0:2048),i_nn(0:2048),j_nn(0:2048),k_nn(0:2048)

!--------------------------------
end module CFD_solver_constants

 !--------------------------------------------------------------------------------------------------------
 !--------------------------------------------------------------------------------------------------------
module  CFD_solver_parameters
   use CFD_solver_constants
   implicit none

! logic------------------------------------------------------------------------------------------
! integer parameters for geometry
   integer,save:: Decomposition, Scheme_type
   integer,save:: nx_global, ny_global, nz_global, nx, ny, nz, LAP
   integer,save:: Dimension, Boundary(3)

!!--------Repeated solver calls
   integer,save:: Istep
   integer,save:: Repeat_count
   real(kind=CFD_solver_real_kind),save:: rparameters(100,10)
   integer,save:: nparameters(100,10)

!-----------------------------------------------------------
! Coordinate parameters
   real(kind=CFD_solver_real_kind),save:: SLx,SLy,SLz,hx,hy,hz

end module  CFD_solver_parameters
!!----------------------------------------------------------


