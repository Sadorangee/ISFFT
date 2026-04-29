module ISFFT_Precision
    implicit none
    include "mpif.h"

    !------For Doubleprecision  (real*8)------------------------
    integer,parameter::ISFFT_REAL_KIND=8,  ISFFT_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
    integer,parameter::ISFFT_DATA_TYPE0=MPI_INTEGER
    ! ------For Single precision (real*4)-----------------------
    ! integer,parameter::ISFFT_REAL_KIND=4,  ISFFT_DATA_TYPE=MPI_REAL             !  single precision computing

 end module ISFFT_Precision