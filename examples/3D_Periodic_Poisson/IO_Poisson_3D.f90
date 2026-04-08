!!------------------------------------------------
subroutine IO_Poisson_3D(phi)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer i,j,k,m, ierr, l
   integer i0,j0,k0
   real(kind=TMS_REAL_KIND)::phi( 1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
   real(kind=TMS_REAL_KIND),allocatable,dimension(:,:,:):: phi_0,phi_00
   character(len=120) :: filename1

   allocate(phi_0( 1:nx_global,1:ny_global,1:nz_global), &
      phi_00(1:nx_global,1:ny_global,1:nz_global))

   phi_0 =0.d0;
   phi_00=0.d0;

!!-------- 3D ----------------


   l= nx_global*ny_global*nz_global;

   DO k= 1,nz
      DO j= 1,ny
         DO i= 1,nx
            i0=i_offset(npx)+i-1;
            j0=j_offset(npy)+j-1;
            k0=k_offset(npz)+k-1;
            phi_0(i0,j0,k0)= phi(i,j,k)
         ENDDO
      ENDDO
   ENDDO

   call MPI_REDUCE(phi_0,phi_00,l,TMS_DATA_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if (my_id == 0) then
      filename1 = 'Poisson_3D_IO.dat'
      print *, 'write data file:', trim(filename1)
      open(153, file=filename1, form='formatted')
      do k = 1, nz_global
         do j = 1, ny_global
            do i = 1, nx_global
               write(153, "(1F16.8)") phi_00(i,j,k)
            end do
         end do
      end do
      close(153)
   end if


end subroutine IO_Poisson_3D
!!------------------------------------------------

