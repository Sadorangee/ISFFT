module  SPLASH_Parameters
   use SPLASH_Precision
   implicit none

   ! ----------------------------------------------------------
   ! integer parameters for geometry
   integer,save:: If_3dfft_decomp, If_scheme
   integer,save:: nx_global,ny_global,nz_global,nx,ny,nz,LAP

   !-----------------------------------------------------------
   ! Coordinate parameters
   real(kind=SPLASH_REAL_KIND),save:: SLx,SLy,SLz

end module  SPLASH_Parameters





