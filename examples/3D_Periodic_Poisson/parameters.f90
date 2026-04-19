!!---------------------------------------------------------------------------------------------------------
module TMS_precision
   implicit none
   include "mpif.h"
!------For Doubleprecision  (real*8)------------------------------------------------------------------
   integer,parameter::TMS_REAL_KIND=8,  TMS_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
   integer,parameter::TMS_DATA_TYPE0=MPI_INTEGER
! ------For Single precision (real*4)-----------------------
   ! integer,parameter::TMS_REAL_KIND=4,  TMS_DATA_TYPE=MPI_REAL             !  single precision computing
!===========  Parameters for MPI ==========================================================----------------
end module TMS_precision

!!------parameters used in TMS solver---------------------------------------
module TMS_constants
!!---------------------------------------------------------------------------------
   Use TMS_precision
   implicit none
!!-------------------------------------------------------------------------
   real(kind=TMS_REAL_KIND),parameter:: TOL_matrix=1e-6;
!--------MPI-------------------
   integer,save:: my_id,npx,npy,npz,npx0,npy0,npz0, ID_XP1,ID_XM1,ID_YP1,ID_YM1,ID_ZP1,ID_ZM1, &
       TMS_Barrier_level,np_size
   integer,save::  i_offset(0:2048),j_offset(0:2048),k_offset(0:2048),i_nn(0:2048),j_nn(0:2048),k_nn(0:2048)

!--------------------------------
end module TMS_constants

 !--------------------------------------------------------------------------------------------------------
 !--------------------------------------------------------------------------------------------------------
module  flow_parameters
   use TMS_constants
   implicit none

! logic------------------------------------------------------------------------------------------
! integer parameters for geometry
   integer,save:: If_3dfft_decomp, If_scheme
   integer,save:: nx_global,ny_global,nz_global,nx,ny,nz,LAP

!!--------Time step
   integer,save:: Istep
   integer,save:: Nt
   real(kind=TMS_REAL_KIND),save:: CFL
   real(kind=TMS_REAL_KIND),save:: dt, T, rparameters(100,10)
   integer,save:: nparameters(100,10)

!-----------------------------------------------------------
! Coordinate parameters
   real(kind=TMS_REAL_KIND),save:: SLx,SLy,SLz,hx,hy,hz
   real(kind=TMS_REAL_KIND),allocatable,save,dimension(:) :: xx,yy,zz,xx0,yy0,zz0

end module  flow_parameters
!!----------------------------------------------------------


