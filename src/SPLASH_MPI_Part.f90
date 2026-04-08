module SPLASH_MPI_Part
   use SPLASH_Parameters
   use SPLASH_MPI_Constants
   implicit none

contains

!----------------------------------------------------------------------
   integer function SPLASH_mod(i,n)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants

      integer, intent(in) :: i,n

      if(i.lt.0) then
         SPLASH_mod=i+n
      else if (i.gt.n-1) then
         SPLASH_mod=i-n
      else
         SPLASH_mod=i
      endif

   end function SPLASH_mod
   !----------------------------------------------------------------------
   subroutine SPLASH_get_i_node(i_global,node_i,i_local)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants

      integer :: i_global,node_i,i_local,ia

      node_i=npx0-1
      do ia=0,npx0-2
         if(i_global.ge.i_offset(ia) .and. i_global .lt. i_offset(ia+1) ) node_i=ia
      enddo
      i_local=i_global-i_offset(node_i)+1

   end subroutine SPLASH_get_i_node
   !----------------------------------------------------------------------
   subroutine SPLASH_get_j_node(j_global,node_j,j_local)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants

      integer :: j_global,node_j,j_local,ja

      node_j=npy0-1
      do ja=0,npy0-2
         if(j_global.ge.j_offset(ja) .and. j_global .lt. j_offset(ja+1) ) node_j=ja
      enddo
      j_local=j_global-j_offset(node_j)+1
   end subroutine SPLASH_get_j_node
   !----------------------------------------------------------------------
   subroutine SPLASH_get_k_node(k_global,node_k,k_local)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants

      integer :: k_global,node_k,k_local,ka

      node_k=npz0-1
      do ka=0,npz0-2
         if(k_global .ge. k_offset(ka) .and. k_global .lt. k_offset(ka+1) ) node_k=ka
      enddo
      k_local=k_global-k_offset(node_k)+1
   end subroutine SPLASH_get_k_node
   !----------------------------------------------------------------------
   integer function SPLASH_get_id(npx1,npy1,npz1)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants

      integer :: npx1,npy1,npz1

      SPLASH_get_id = npz1*(npx0*npy0)+npy1*npx0+npx1

   end function SPLASH_get_id
   !----------------------------------------------------------------------


   !----------------------------------------------------------------------
   subroutine SPLASH_part()
      use SPLASH_Parameters
      use SPLASH_MPI_Constants
      
      integer :: k,ka,ierr
      integer :: npx1,npy1,npz1,npx2,npy2,npz2

      if(np_size .ne. npx0*npy0*npz0) then
         if(my_id.eq.0) print*, 'The Number of total Processes is not equal to npx0*npy0*npz0 !'
         call mpi_finalize(ierr)
         stop
      endif

      npx=mod(my_id,npx0)
      npy=mod(my_id,npx0*npy0)/npx0
      npz=my_id/(npx0*npy0)

      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npy*npx0, npx,MPI_COMM_X,ierr)   ! 1-D
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npx, npy,MPI_COMM_Y,ierr)
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npy*npx0+npx, npz,MPI_COMM_Z,ierr)

      nx=nx_global/npx0
      ny=ny_global/npy0
      nz=nz_global/npz0

      if(npx .lt. mod(nx_global,npx0)) nx=nx+1
      if(npy .lt. mod(ny_global,npy0)) ny=ny+1
      if(npz .lt. mod(nz_global,npz0)) nz=nz+1

      do k=0,npx0-1
         ka=min(k,mod(nx_global,npx0))
         i_offset(k)=int(nx_global/npx0)*k+ka+1
         i_nn(k)=nx_global/npx0
         if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
      enddo

      do k=0,npy0-1
         ka=min(k,mod(ny_global,npy0))
         j_offset(k)=int(ny_global/npy0)*k+ka+1
         j_nn(k)=ny_global/npy0
         if(k .lt. mod(ny_global,npy0)) j_nn(k)=j_nn(k)+1
      enddo

      do k=0,npz0-1
         ka=min(k,mod(nz_global,npz0))
         k_offset(k)=int(nz_global/npz0)*k+ka+1
         k_nn(k)=nz_global/npz0
         if(k .lt. mod(nz_global,npz0)) k_nn(k)=k_nn(k)+1
      enddo

      npx1=SPLASH_mod(npx-1,npx0)
      npx2=SPLASH_mod(npx+1,npx0)
      ID_XM1=npz*(npx0*npy0)+npy*npx0+npx1    ! -1 proc in x-direction
      ID_XP1=npz*(npx0*npy0)+npy*npx0+npx2    ! +1 proc in x-direction


      npy1=SPLASH_mod(npy-1,npy0)
      npy2=SPLASH_mod(npy+1,npy0)
      ID_YM1=npz*(npx0*npy0)+npy1*npx0+npx
      ID_YP1=npz*(npx0*npy0)+npy2*npx0+npx


      npz1=SPLASH_mod(npz-1,npz0)
      npz2=SPLASH_mod(npz+1,npz0)
      ID_ZM1=npz1*(npx0*npy0)+npy*npx0+npx
      ID_ZP1=npz2*(npx0*npy0)+npy*npx0+npx

      call MPI_barrier(MPI_COMM_WORLD,ierr)

   end subroutine SPLASH_part
   !----------------------------------------------------------------------
   subroutine SPLASH_part_change(npx0_new,npy0_new,npz0_new)
      use SPLASH_Parameters
      use SPLASH_MPI_Constants
      integer :: npx0_new, npy0_new, npz0_new
      integer :: ierr

      npx0 = npx0_new;
      npy0 = npy0_new;
      npz0 = npz0_new;

      if(np_size .ne. npx0*npy0*npz0) then
         if(my_id.eq.0) print*, 'The Changed Number of total Processes is not equal to npx0*npy0*npz0 !'
         call mpi_finalize(ierr)
         stop
      endif

      call mpi_comm_free(MPI_COMM_X,ierr)
      call mpi_comm_free(MPI_COMM_Y,ierr)
      call mpi_comm_free(MPI_COMM_Z,ierr)
      call SPLASH_part

   end subroutine SPLASH_part_change
   !----------------------------------------------------------------------


end module SPLASH_MPI_Part
