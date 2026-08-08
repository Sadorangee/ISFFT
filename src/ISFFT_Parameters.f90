module  ISFFT_Parameters
   use ISFFT_Precision
   implicit none

   ! ----------------------------------------------------------
   ! integer parameters for geometry
   integer,save:: If_decomp, If_scheme, Dimension
   integer,save:: ISFFT_plan_generation = 0
   integer,save:: nx_global,ny_global,nz_global,nx,ny,nz,LAP
   integer,save:: BC(3)

   !-----------------------------------------------------------
   ! Coordinate parameters
   real(kind=ISFFT_REAL_KIND),save:: SLx,SLy,SLz

end module  ISFFT_Parameters


