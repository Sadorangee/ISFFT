module SPLASH_Precision
    implicit none
    include "mpif.h"

    !------For Doubleprecision  (real*8)------------------------
    integer,parameter::SPLASH_REAL_KIND=8,  SPLASH_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
    integer,parameter::SPLASH_DATA_TYPE0=MPI_INTEGER
    ! ------For Single precision (real*4)-----------------------
    ! integer,parameter::SPLASH_REAL_KIND=4,  SPLASH_DATA_TYPE=MPI_REAL             !  single precision computing

 end module SPLASH_Precision