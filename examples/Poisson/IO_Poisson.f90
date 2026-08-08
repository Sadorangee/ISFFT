!!------------------------------------------------
subroutine IO_Poisson_1D(phi)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer :: i
   real(kind=TMS_REAL_KIND) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)

   if (my_id == 0) then
      print *, 'write data file: Poisson_1D_IO.dat'
      open(153, file='Poisson_1D_IO.dat', form='formatted')
      do i = 1, nx_global
         write(153, '(1F24.16)') phi(i,1,1)
      end do
      close(153)
   end if
end subroutine IO_Poisson_1D
!!------------------------------------------------

!!------------------------------------------------
subroutine IO_Poisson_2D(phi)
   Use flow_parameters
   use TMS_constants
   implicit none

   integer :: i, j, i0, j0, ierr
   real(kind=TMS_REAL_KIND) :: phi(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
   real(kind=TMS_REAL_KIND), allocatable :: phi_local(:,:), phi_global(:,:)

   allocate(phi_local(nx_global,ny_global), source=0.0_TMS_REAL_KIND)
   allocate(phi_global(nx_global,ny_global), source=0.0_TMS_REAL_KIND)
   do j = 1, ny
      j0 = j_offset(npy)+j-1
      do i = 1, nx
         i0 = i_offset(npx)+i-1
         phi_local(i0,j0) = phi(i,j,1)
      end do
   end do
   call MPI_Reduce(phi_local, phi_global, nx_global*ny_global, TMS_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

   if (my_id == 0) then
      print *, 'write data file: Poisson_2D_IO.dat'
      open(153, file='Poisson_2D_IO.dat', form='formatted')
      do j = 1, ny_global
         do i = 1, nx_global
            write(153, '(1F24.16)') phi_global(i,j)
         end do
      end do
      close(153)
   end if
   deallocate(phi_local, phi_global)
end subroutine IO_Poisson_2D
!!------------------------------------------------

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

!!!------------------------------------------------

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
               write(153, '(1F24.16)') phi_00(i,j,k)
            end do
         end do
      end do
      close(153)
   end if


end subroutine IO_Poisson_3D
!!------------------------------------------------

