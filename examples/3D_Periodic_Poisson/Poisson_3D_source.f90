!!-----------------------------------------
subroutine Poisson_3D_source(f)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer i,j,k,m, ierr, l, k1
   real(kind=TMS_REAL_KIND)::f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

   real(kind=TMS_REAL_KIND),parameter:: pi=3.14159260d0;
   real(kind=TMS_REAL_KIND) :: sx, sy, sz

   allocate(xx(1-LAP:nx+LAP),yy(1-LAP:ny+LAP),zz(1-LAP:nz+LAP),  &
      xx0(nx_global), yy0(ny_global), zz0(nz_global) )

   f=0.0;
   xx =0.0;
   yy =0.0;
   zz =0.0;
   xx0=0.0;
   yy0=0.0;
   zz0=0.0;

   sx = pi/SLx;
   sy = pi/SLy;
   sz = pi/SLz;

   do i=1,nx_global
      xx0(i)=(i-1.0)*hx
   enddo
   do k=1-LAP,nx+LAP
      k1=i_offset(npx)+k-1
      if(k1.ge.1 .and. k1 .le. nx_global) then
         xx(k)=xx0(k1)
      endif
   enddo
   do j=1,ny_global
      yy0(j)=(j-1.0)*hy
   enddo
   do k=1-LAP,ny+LAP
      k1=j_offset(npy)+k-1
      if(k1.ge.1 .and. k1 .le. ny_global) then
         yy(k)=yy0(k1)
      endif
   enddo
   do k=1,nz_global
      zz0(k)=(k-1.0)*hz
   enddo
   do k=1-LAP,nz+LAP
      k1=k_offset(npz)+k-1
      if(k1.ge.1 .and. k1 .le. nz_global) then
         zz(k)=zz0(k1)
      endif
   enddo

   !! 3D Periodic Poisson 
   do k= 1,nz
      do j= 1,ny
         do i= 1,nx
            f(i,j,k)= -sin(2.0*sz*zz(k))*(((2.0*sx)**2.0+(2.0*sz)**2.0)*sin(2.0*sx*xx(i)) + &
                      ((4.0*sy)**2.0+(2.0*sz)**2.0)*cos(4.0*sy*yy(j)))
         enddo
      enddo
   enddo

   deallocate(xx,yy,zz,xx0,yy0,zz0)

end  subroutine Poisson_3D_source
!!-----------------------------------------

