!!-----------------------------------------
subroutine Poisson_source(f)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer i, j, k, m, ierr, l, k1
   real(kind=TMS_REAL_KIND)::f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

   real(kind=TMS_REAL_KIND),parameter:: PI=3.14159265358979323846d0;
   real(kind=TMS_REAL_KIND) :: sx, sy, sz
   real(kind=TMS_REAL_KIND) :: gx, gy, gz, g2x, g2y,g2z
   real(kind=TMS_REAL_KIND) :: phix, phiy, phiz, mux, muy, muz

   real(kind=TMS_REAL_KIND) :: sigma
   real(kind=TMS_REAL_KIND) :: xc, yc, zc
   real(kind=TMS_REAL_KIND) :: dx, dy, dz
   real(kind=TMS_REAL_KIND) :: ex, ey, ez


   real(kind=TMS_REAL_KIND) :: s_width
   real(kind=TMS_REAL_KIND) :: xstart, xend, ystart, yend, zstart, zend
   real(kind=TMS_REAL_KIND) :: ux, vx, uy, vy, uz, vz
   real(kind=TMS_REAL_KIND) :: txu, txv, tyu, tyv, tzu, tzv
   real(kind=TMS_REAL_KIND) :: sxu2, sxv2, syu2, syv2, szu2, szv2

   integer :: isrc, jsrc, ksrc
   integer :: ig, jg, kg
   real(kind=TMS_REAL_KIND) :: dV


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

   if (Dimension == 1) then
      do i = 1, nx
         select case (Boundary(1))
         case (0)
            ! 1D P case 0: phi = sin(2*pi*x/Lx)
            f(i,1,1) = -(2.0d0*sx)**2 * sin(2.0d0*sx*xx(i))
         case (1)
            ! 1D D case 0: phi = sin(pi*x/Lx)
            f(i,1,1) = -sx**2 * sin(sx*xx(i))
         case (2)
            ! 1D N case 0: phi = cos(pi*x/Lx)
            f(i,1,1) = -sx**2 * cos(sx*xx(i))
         end select
      end do
      deallocate(xx,yy,zz,xx0,yy0,zz0)
      return
   end if

   if (Dimension == 2) then
      do j = 1, ny
         do i = 1, nx
            select case (10*Boundary(1)+Boundary(2))
            case (0)
               ! 2D PP case 0: phi = sin(2*pi*x/Lx) + cos(4*pi*y/Ly)
               f(i,j,1) = -(2.0d0*sx)**2*sin(2.0d0*sx*xx(i)) &
                  -(4.0d0*sy)**2*cos(4.0d0*sy*yy(j))
            case (11)
               ! 2D DD case 0: phi = sin(pi*x/Lx) sin(2*pi*y/Ly)
               f(i,j,1) = -(sx**2+(2.0d0*sy)**2) &
                  *sin(sx*xx(i))*sin(2.0d0*sy*yy(j))
            case (22)
               ! 2D NN case 0: phi = cos(pi*x/Lx) cos(2*pi*y/Ly)
               f(i,j,1) = -(sx**2+(2.0d0*sy)**2) &
                  *cos(sx*xx(i))*cos(2.0d0*sy*yy(j))
            case (1)
               ! 2D PD case 0: phi = cos(2*pi*x/Lx) sin(pi*y/Ly)
               f(i,j,1) = -((2.0d0*sx)**2+sy**2) &
                  *cos(2.0d0*sx*xx(i))*sin(sy*yy(j))
            case (2)
               ! 2D PN case 0: phi = cos(2*pi*x/Lx) cos(pi*y/Ly)
               f(i,j,1) = -((2.0d0*sx)**2+sy**2) &
                  *cos(2.0d0*sx*xx(i))*cos(sy*yy(j))
            case (12)
               ! 2D DN case 0: phi = sin(pi*x/Lx) cos(pi*y/Ly)
               f(i,j,1) = -(sx**2+sy**2)*sin(sx*xx(i))*cos(sy*yy(j))
            case default
               ! Reversed mixed order remains available for solver regression cases.
               if (Boundary(1) == 0) then
                  phix = cos(2.0d0*sx*xx(i)); mux = (2.0d0*sx)**2
               else if (Boundary(1) == 1) then
                  phix = sin(sx*xx(i)); mux = sx**2
               else
                  phix = cos(sx*xx(i)); mux = sx**2
               end if
               if (Boundary(2) == 0) then
                  phiy = cos(2.0d0*sy*yy(j)); muy = (2.0d0*sy)**2
               else if (Boundary(2) == 1) then
                  phiy = sin(sy*yy(j)); muy = sy**2
               else
                  phiy = cos(sy*yy(j)); muy = sy**2
               end if
               f(i,j,1) = -(mux+muy)*phix*phiy
            end select
         end do
      end do
      deallocate(xx,yy,zz,xx0,yy0,zz0)
      return
   end if

   if (If_Dirichlet_test == 7) then
         dV = hx * hy * hz

         isrc = nx_global / 2
         jsrc = ny_global / 2
         ksrc = nz_global / 2

         isrc = max(2, min(nx_global-1, isrc))
         jsrc = max(2, min(ny_global-1, jsrc))
         ksrc = max(2, min(nz_global-1, ksrc))
   endif

   !! 3D Periodic Poisson
   if ((Boundary(1) .eq. 0) .and. (Boundary(2) .eq. 0) .and. (Boundary(3) .eq. 0)) then
      if (If_Periodic_test == 0) then
         do k= 1,nz
            do j= 1,ny
               do i= 1,nx
                  f(i,j,k)= -sin(2.0*sz*zz(k))*(((2.0*sx)**2.0+(2.0*sz)**2.0)*sin(2.0*sx*xx(i)) + &
                     ((4.0*sy)**2.0+(2.0*sz)**2.0)*cos(4.0*sy*yy(j)))
               enddo
            enddo
         enddo
      endif

      if (If_Periodic_test == 1) then
         do k= 1,nz
            do j= 1,ny
               do i= 1,nx
                  !! Polynomial source for periodic BCs
                  !! phi(x,y,z) = x^2(SLx-x)^2 * y^2(SLy-y)^2 * z^2(SLz-z)^2
                  gx = (xx(i)**2) * ((SLx-xx(i))**2)
                  gy = (yy(j)**2) * ((SLy-yy(j))**2)
                  gz = (zz(k)**2) * ((SLz-zz(k))**2)

                  g2x = 12.0d0*(xx(i)**2) - 12.0d0*SLx*xx(i) + 2.0d0*(SLx**2)
                  g2y = 12.0d0*(yy(j)**2) - 12.0d0*SLy*yy(j) + 2.0d0*(SLy**2)
                  g2z = 12.0d0*(zz(k)**2) - 12.0d0*SLz*zz(k) + 2.0d0*(SLz**2)

                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               enddo
            enddo
         enddo
      endif
   end if

   !! 3D mixed-boundary manufactured Poisson case.
   !! Periodic  : cos(2*pi*x/L)
   !! Dirichlet : sin(pi*x/L)
   !! Neumann   : cos(pi*x/L)
   if (.not. ((Boundary(1) .eq. Boundary(2)) .and. (Boundary(2) .eq. Boundary(3)))) then
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (Boundary(1) .eq. 0) then
                  phix = cos(2.0d0*pi*xx(i)/SLx)
                  mux = (2.0d0*pi/SLx)**2
               else if (Boundary(1) .eq. 1) then
                  phix = sin(pi*xx(i)/SLx)
                  mux = (pi/SLx)**2
               else
                  phix = cos(pi*xx(i)/SLx)
                  mux = (pi/SLx)**2
               end if

               if (Boundary(2) .eq. 0) then
                  phiy = cos(2.0d0*pi*yy(j)/SLy)
                  muy = (2.0d0*pi/SLy)**2
               else if (Boundary(2) .eq. 1) then
                  phiy = sin(pi*yy(j)/SLy)
                  muy = (pi/SLy)**2
               else
                  phiy = cos(pi*yy(j)/SLy)
                  muy = (pi/SLy)**2
               end if

               if (Boundary(3) .eq. 0) then
                  phiz = cos(2.0d0*pi*zz(k)/SLz)
                  muz = (2.0d0*pi/SLz)**2
               else if (Boundary(3) .eq. 1) then
                  phiz = sin(pi*zz(k)/SLz)
                  muz = (pi/SLz)**2
               else
                  phiz = cos(pi*zz(k)/SLz)
                  muz = (pi/SLz)**2
               end if

               f(i,j,k) = -(mux + muy + muz)*phix*phiy*phiz
            end do
         end do
      end do
   end if

   !! 3D Dirichlet Poisson
   if ((Boundary(1) .eq. 1) .and. (Boundary(2) .eq. 1) .and. (Boundary(3) .eq. 1)) then
      do k= 1,nz
         do j= 1,ny
            do i= 1,nx
               if (If_Dirichlet_test == 0) f(i,j,k) = - ( sx**2 + (2.0d0*sy)**2 + (3.0d0*sz)**2 ) * &
                  sin(sx*xx(i)) * sin(2.0d0*sy*yy(j)) * sin(3.0d0*sz*zz(k))
               if (If_Dirichlet_test == 1) f(i,j,k) = -3.0d0 * pi **2.0 * &
                  sin(pi*xx(i))*sin(pi*yy(j))*sin(pi*zz(k))
               if (If_Dirichlet_test == 2) f(i,j,k) = -2.0d0 * ( &
                  (yy(j)*(SLy-yy(j))) * (zz(k)*(SLz-zz(k))) + &
                  (xx(i)*(SLx-xx(i))) * (zz(k)*(SLz-zz(k))) + &
                  (xx(i)*(SLx-xx(i))) * (yy(j)*(SLy-yy(j))) )
               if (If_Dirichlet_test == 3) then
                  ! Exact solution:
                  ! phi = [x^2 (Lx-x)^2] [y^2 (Ly-y)^2] [z^2 (Lz-z)^2]

                  gx  = (xx(i)**2) * ((SLx-xx(i))**2)
                  gy  = (yy(j)**2) * ((SLy-yy(j))**2)
                  gz  = (zz(k)**2) * ((SLz-zz(k))**2)

                  g2x = 12.0d0*(xx(i)**2) - 12.0d0*SLx*xx(i) + 2.0d0*(SLx**2)
                  g2y = 12.0d0*(yy(j)**2) - 12.0d0*SLy*yy(j) + 2.0d0*(SLy**2)
                  g2z = 12.0d0*(zz(k)**2) - 12.0d0*SLz*zz(k) + 2.0d0*(SLz**2)

                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               endif
               if (If_Dirichlet_test == 4) then
                  ! Exact solution:
                  ! phi = [x^3 (Lx-x)^3] [y^3 (Ly-y)^3] [z^3 (Lz-z)^3]

                  gx  = (xx(i)**3) * ((SLx-xx(i))**3)
                  gy  = (yy(j)**3) * ((SLy-yy(j))**3)
                  gz  = (zz(k)**3) * ((SLz-zz(k))**3)

                  g2x = 6.0d0*(SLx**3)*xx(i) - 36.0d0*(SLx**2)*(xx(i)**2) + &
                     60.0d0*SLx*(xx(i)**3) - 30.0d0*(xx(i)**4)

                  g2y = 6.0d0*(SLy**3)*yy(j) - 36.0d0*(SLy**2)*(yy(j)**2) + &
                     60.0d0*SLy*(yy(j)**3) - 30.0d0*(yy(j)**4)

                  g2z = 6.0d0*(SLz**3)*zz(k) - 36.0d0*(SLz**2)*(zz(k)**2) + &
                     60.0d0*SLz*(zz(k)**3) - 30.0d0*(zz(k)**4)

                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               endif
               if (If_Dirichlet_test == 5) then
                  ! Exact solution:
                  ! phi(x,y,z) =
                  ! exp(-((x-SLx/2)^2)/(2*sigma^2)) * sin(pi*x/SLx)^2
                  ! * exp(-((y-SLy/2)^2)/(2*sigma^2)) * sin(pi*y/SLy)^2
                  ! * exp(-((z-SLz/2)^2)/(2*sigma^2)) * sin(pi*z/SLz)^2
                  !
                  ! Poisson equation:
                  ! laplacian(phi) = f
                  !
                  ! where
                  ! f = g2x*gy*gz + gx*g2y*gz + gx*gy*g2z

                  sigma = 0.4d0

                  xc = 0.5d0 * SLx
                  yc = 0.5d0 * SLy
                  zc = 0.5d0 * SLz

                  dx = xx(i) - xc
                  dy = yy(j) - yc
                  dz = zz(k) - zc

                  ex = exp( - dx**2 / (2.0d0*sigma**2) )
                  ey = exp( - dy**2 / (2.0d0*sigma**2) )
                  ez = exp( - dz**2 / (2.0d0*sigma**2) )

                  gx = ex * ( sin(pi*xx(i)/SLx) )**2
                  gy = ey * ( sin(pi*yy(j)/SLy) )**2
                  gz = ez * ( sin(pi*zz(k)/SLz) )**2

                  g2x = ex * ( &
                     ( dx**2 / sigma**4 - 1.0d0 / sigma**2 ) * ( sin(pi*xx(i)/SLx) )**2 &
                     - ( 2.0d0*pi*dx / ( sigma**2 * SLx ) ) * sin( 2.0d0*pi*xx(i)/SLx ) &
                     + 2.0d0 * ( pi / SLx )**2 * cos( 2.0d0*pi*xx(i)/SLx ) )

                  g2y = ey * ( &
                     ( dy**2 / sigma**4 - 1.0d0 / sigma**2 ) * ( sin(pi*yy(j)/SLy) )**2 &
                     - ( 2.0d0*pi*dy / ( sigma**2 * SLy ) ) * sin( 2.0d0*pi*yy(j)/SLy ) &
                     + 2.0d0 * ( pi / SLy )**2 * cos( 2.0d0*pi*yy(j)/SLy ) )

                  g2z = ez * ( &
                     ( dz**2 / sigma**4 - 1.0d0 / sigma**2 ) * ( sin(pi*zz(k)/SLz) )**2 &
                     - ( 2.0d0*pi*dz / ( sigma**2 * SLz ) ) * sin( 2.0d0*pi*zz(k)/SLz ) &
                     + 2.0d0 * ( pi / SLz )**2 * cos( 2.0d0*pi*zz(k)/SLz ) )

                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               endif
               if (If_Dirichlet_test == 6) then
                  ! Exact solution:
                  ! phi(x,y,z) = gx(x) * gy(y) * gz(z)
                  !
                  ! gx(x) = 0.5*(1 + tanh((x-xstart)/s)) * 0.5*(1 - tanh((x-xend)/s))
                  ! gy(y) = 0.5*(1 + tanh((y-ystart)/s)) * 0.5*(1 - tanh((y-yend)/s))
                  ! gz(z) = 0.5*(1 + tanh((z-zstart)/s)) * 0.5*(1 - tanh((z-zend)/s))
                  !
                  ! Poisson equation:
                  ! laplacian(phi) = f
                  !
                  ! Note:
                  ! this function is near-zero at boundaries, but not exactly zero in general.
               
                  s_width = 0.2d0
               
                  xstart = 0.25d0 * SLx
                  xend   = 0.75d0 * SLx
                  ystart = 0.25d0 * SLy
                  yend   = 0.75d0 * SLy
                  zstart = 0.25d0 * SLz
                  zend   = 0.75d0 * SLz
               
                  ux = (xx(i) - xstart) / s_width
                  vx = (xx(i) - xend  ) / s_width
                  uy = (yy(j) - ystart) / s_width
                  vy = (yy(j) - yend  ) / s_width
                  uz = (zz(k) - zstart) / s_width
                  vz = (zz(k) - zend  ) / s_width
               
                  txu = tanh(ux)
                  txv = tanh(vx)
                  tyu = tanh(uy)
                  tyv = tanh(vy)
                  tzu = tanh(uz)
                  tzv = tanh(vz)
               
                  sxu2 = 1.0d0 / cosh(ux)**2
                  sxv2 = 1.0d0 / cosh(vx)**2
                  syu2 = 1.0d0 / cosh(uy)**2
                  syv2 = 1.0d0 / cosh(vy)**2
                  szu2 = 1.0d0 / cosh(uz)**2
                  szv2 = 1.0d0 / cosh(vz)**2
               
                  gx = 0.25d0 * (1.0d0 + txu) * (1.0d0 - txv)
                  gy = 0.25d0 * (1.0d0 + tyu) * (1.0d0 - tyv)
                  gz = 0.25d0 * (1.0d0 + tzu) * (1.0d0 - tzv)
               
                  g2x = -(0.5d0 / s_width**2) * (1.0d0 - txv) * sxu2 * txu &
                        -(0.5d0 / s_width**2) * sxu2 * sxv2 &
                        +(0.5d0 / s_width**2) * (1.0d0 + txu) * sxv2 * txv
               
                  g2y = -(0.5d0 / s_width**2) * (1.0d0 - tyv) * syu2 * tyu &
                        -(0.5d0 / s_width**2) * syu2 * syv2 &
                        +(0.5d0 / s_width**2) * (1.0d0 + tyu) * syv2 * tyv
               
                  g2z = -(0.5d0 / s_width**2) * (1.0d0 - tzv) * szu2 * tzu &
                        -(0.5d0 / s_width**2) * szu2 * szv2 &
                        +(0.5d0 / s_width**2) * (1.0d0 + tzu) * szv2 * tzv
               
                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               endif

               if (If_Dirichlet_test == 7) then
                  ig = i_offset(npx) + i - 1
                  jg = j_offset(npy) + j - 1
                  kg = k_offset(npz) + k - 1

                  if ((ig == isrc) .and. (jg == jsrc) .and. (kg == ksrc)) then
                     f(i,j,k) = 1.0d0 / dV
                  else
                     f(i,j,k) = 0.0d0
                  endif
               endif

            enddo
         enddo
      enddo
   end if

   if ((Boundary(1) .eq. 2) .and. (Boundary(2) .eq. 2) .and. (Boundary(3) .eq. 2)) then

      if (If_Neumann_test == 2) then
         dV = hx * hy * hz

         isrc = nx_global / 2
         jsrc = ny_global / 2
         ksrc = nz_global / 2

         isrc = max(1, min(nx_global, isrc))
         jsrc = max(1, min(ny_global, jsrc))
         ksrc = max(1, min(nz_global, ksrc))
      endif

      do k= 1,nz
         do j= 1,ny
            do i= 1,nx

               if (If_Neumann_test == 0) then
                  ! ---------------------------------------------------------
                  ! Case 0: Trigonometric Neumann case
                  !
                  ! Exact solution:
                  ! phi = cos(pi*x) cos(2*pi*y) cos(3*pi*z)
                  !
                  ! Neumann BC:
                  ! dphi/dn = 0 on all boundaries
                  !
                  ! Source:
                  ! laplacian(phi) = f
                  ! ---------------------------------------------------------
                  f(i,j,k) = - ( pi**2 + (2.0d0*pi)**2 + (3.0d0*pi)**2 ) * &
                                 cos(pi*xx(i)) * cos(2.0d0*pi*yy(j)) * cos(3.0d0*pi*zz(k))
               endif


               if (If_Neumann_test == 1) then
                  ! ---------------------------------------------------------
                  ! Case 1: Polynomial Neumann case
                  !
                  ! Exact solution:
                  ! phi = gx(x) gy(y) gz(z)
                  !
                  ! gx = x^2 (Lx-x)^2
                  ! gy = y^2 (Ly-y)^2
                  ! gz = z^2 (Lz-z)^2
                  !
                  ! This satisfies:
                  ! dphi/dn = 0 on all boundaries
                  !
                  ! Source:
                  ! f = g2x*gy*gz + gx*g2y*gz + gx*gy*g2z
                  ! ---------------------------------------------------------

                  gx  = (xx(i)**2) * ((SLx-xx(i))**2)
                  gy  = (yy(j)**2) * ((SLy-yy(j))**2)
                  gz  = (zz(k)**2) * ((SLz-zz(k))**2)

                  g2x = 12.0d0*(xx(i)**2) - 12.0d0*SLx*xx(i) + 2.0d0*(SLx**2)
                  g2y = 12.0d0*(yy(j)**2) - 12.0d0*SLy*yy(j) + 2.0d0*(SLy**2)
                  g2z = 12.0d0*(zz(k)**2) - 12.0d0*SLz*zz(k) + 2.0d0*(SLz**2)

                  f(i,j,k) = g2x * gy * gz + gx * g2y * gz + gx * gy * g2z
               endif


               if (If_Neumann_test == 2) then
                  ! ---------------------------------------------------------
                  ! Case 2: Discrete Dirac source for pure Neumann Poisson
                  !
                  ! Continuous form:
                  ! laplacian(phi) = delta(x-x0,y-y0,z-z0) - 1/|Omega|
                  !
                  ! For pure Neumann, the source must have zero mean.
                  ! Here a discrete zero-mean version is used:
                  !
                  ! f = delta_h - mean(delta_h)
                  !
                  ! where delta_h = 1/dV at the source grid point.
                  !
                  ! This case has a singular Green-function solution.
                  ! It is not recommended for MAE comparison against a smooth
                  ! exact solution.
                  ! ---------------------------------------------------------

                  ig = i_offset(npx) + i - 1
                  jg = j_offset(npy) + j - 1
                  kg = k_offset(npz) + k - 1

                  f(i,j,k) = - 1.0d0 / (dV * real(nx_global*ny_global*nz_global, TMS_REAL_KIND))

                  if ((ig == isrc) .and. (jg == jsrc) .and. (kg == ksrc)) then
                     f(i,j,k) = f(i,j,k) + 1.0d0 / dV
                  endif
               endif

            enddo
         enddo
      enddo
   end if


   deallocate(xx,yy,zz,xx0,yy0,zz0)

end  subroutine Poisson_source
!!-----------------------------------------
