module  ISFFT_Parameters
   use ISFFT_Precision
   implicit none

   ! ----------------------------------------------------------
   ! integer parameters for geometry
   integer,save:: If_3dfft_decomp, If_scheme
   integer,save:: nx_global,ny_global,nz_global,nx,ny,nz,LAP

   !-----------------------------------------------------------
   ! Coordinate parameters
   real(kind=ISFFT_REAL_KIND),save:: SLx,SLy,SLz

end module  ISFFT_Parameters





