!------------------------------------------------------------------------------
! ISFFT
! Public library imports and exported API.
!------------------------------------------------------------------------------
module ISFFT
   !! Imported public interfaces and constants
   use ISFFT_Problem_State, only: ISFFT_PLAN_ESTIMATE, ISFFT_PLAN_MEASURE, ISFFT_USE_REDISTRIBUTION, ISFFT_USE_POISSON, &
      ISFFT_USE_REDISTRIBUTION_POISSON, ISFFT_USE_ALL, ISFFT_USE_R2CFFT_OUTOFPLACE, ISFFT_USE_R2CFFT_INPLACE, &
      ISFFT_USE_C2CFFT_OUTOFPLACE, ISFFT_USE_C2CFFT_INPLACE, ISFFT_USE_REDISTRIBUTION_R2CFFT_OUTOFPLACE, &
      ISFFT_USE_REDISTRIBUTION_R2CFFT_INPLACE, ISFFT_USE_REDISTRIBUTION_C2CFFT_OUTOFPLACE, &
      ISFFT_USE_REDISTRIBUTION_C2CFFT_INPLACE, ISFFT_FFT_INTERNAL, ISFFT_FFT_ORIGINAL, ISFFT_FFT_BOTH
   use ISFFT_Initialization, only: ISFFT_Set_Plan_Options => Set_Plan_Options
   use ISFFT_Precision, only: ISFFT_REAL_KIND, ISFFT_REAL_MPI_TYPE
   use ISFFT_Initialization, only: ISFFT_Initialize, ISFFT_Finalize, ISFFT_Get_Local_Layout
   use ISFFT_Integrated_Solver, only: ISFFT_Poisson_solver
   use ISFFT_Transform_Execution, only: ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c, ISFFT_2D_FFT_r2c, ISFFT_2D_FFT_c2c, &
      ISFFT_3D_FFT_r2c, ISFFT_3D_FFT_c2c, ISFFT_1D_FFT_r2c_inplace, ISFFT_1D_FFT_c2c_inplace, &
      ISFFT_2D_FFT_r2c_inplace, ISFFT_2D_FFT_c2c_inplace, ISFFT_3D_FFT_r2c_inplace, ISFFT_3D_FFT_c2c_inplace, &
      ISFFT_Get_FFT_Layout, ISFFT_Get_FFT_Local_Sizes
   use ISFFT_Redistribution, only: ISFFT_Redistribution_Forward, ISFFT_Redistribution_Inverse
   implicit none

   private
   !! Exported API
   !! Keep exports aligned with the imports above when editing the public API.
   public :: ISFFT_PLAN_ESTIMATE, ISFFT_PLAN_MEASURE, ISFFT_USE_REDISTRIBUTION, ISFFT_USE_POISSON, &
      ISFFT_USE_REDISTRIBUTION_POISSON, ISFFT_USE_ALL, ISFFT_USE_R2CFFT_OUTOFPLACE, ISFFT_USE_R2CFFT_INPLACE, &
      ISFFT_USE_C2CFFT_OUTOFPLACE, ISFFT_USE_C2CFFT_INPLACE, ISFFT_USE_REDISTRIBUTION_R2CFFT_OUTOFPLACE, &
      ISFFT_USE_REDISTRIBUTION_R2CFFT_INPLACE, ISFFT_USE_REDISTRIBUTION_C2CFFT_OUTOFPLACE, &
      ISFFT_USE_REDISTRIBUTION_C2CFFT_INPLACE, ISFFT_FFT_INTERNAL, ISFFT_FFT_ORIGINAL, ISFFT_FFT_BOTH
   public :: ISFFT_Set_Plan_Options
   public :: ISFFT_REAL_KIND, ISFFT_REAL_MPI_TYPE
   public :: ISFFT_Initialize, ISFFT_Finalize, ISFFT_Get_Local_Layout
   public :: ISFFT_Poisson_solver
   public :: ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c
   public :: ISFFT_2D_FFT_r2c, ISFFT_2D_FFT_c2c
   public :: ISFFT_3D_FFT_r2c, ISFFT_3D_FFT_c2c
   public :: ISFFT_1D_FFT_r2c_inplace, ISFFT_1D_FFT_c2c_inplace
   public :: ISFFT_2D_FFT_r2c_inplace, ISFFT_2D_FFT_c2c_inplace
   public :: ISFFT_3D_FFT_r2c_inplace, ISFFT_3D_FFT_c2c_inplace
   public :: ISFFT_Get_FFT_Layout, ISFFT_Get_FFT_Local_Sizes
   public :: ISFFT_Redistribution_Forward, ISFFT_Redistribution_Inverse

end module ISFFT
