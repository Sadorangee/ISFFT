!!-----------------------------------------
subroutine Poisson_3D_source(f)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer i,j,k,m, ierr, l, k1
   real(kind=TMS_REAL_KIND)::f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

   real(kind=TMS_REAL_KIND),parameter:: PI=3.14159260d0;

   allocate(xx(1-LAP:nx+LAP),yy(1-LAP:ny+LAP),zz(1-LAP:nz+LAP),  &
      xx0(nx_global), yy0(ny_global), zz0(nz_global) )

   allocate(sz(1-LAP:nz+LAP),sz0(nz_global));

   f=0.0;
   xx =0.0;yy =0.0;zz =0.0;
   xx0=0.0;yy0=0.0;zz0=0.0;
   sz =0.d0;
   sz0=0.d0;

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
      sz0(k)=1.0;
   enddo
   do k=1-LAP,nz+LAP
      k1=k_offset(npz)+k-1
      if(k1.ge.1 .and. k1 .le. nz_global) then
         zz(k)=zz0(k1)
         sz(k)=sz0(k1)
      endif
   enddo

   !! 3D Periodic Poisson 
   do k= 1,nz
      do j= 1,ny
         do i= 1,nx
            f(i,j,k)= -pi**2.0*sin(pi*zz(k))*(2.0*sin(pi*xx(i))+5.0*cos(2.0*pi*yy(j)));
         enddo
      enddo
   enddo

         



   deallocate(xx,yy,zz,xx0,yy0,zz0,sz,sz0)

end  subroutine Poisson_3D_source
!!-----------------------------------------

