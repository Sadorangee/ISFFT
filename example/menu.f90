!!--------------------------------------------------------------------------------
include'parameters.f90'
include'part_mpi.f90'
include'main_solver.f90'

Program menu
use CFD_solver_parameters
use CFD_solver_constants
implicit none
integer :: ierr
!---------------------------------------------------------------------------------
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
call mpi_comm_size(MPI_COMM_WORLD,np_size,ierr)
call read_parameters

!------ npx0,npy0,npz0 are parallel partitions in x, y and z directions----------
call part

call main_solver

call mpi_finalize(ierr)

end Program menu
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine read_parameters
   use CFD_solver_parameters
   implicit none
   integer :: ierr, input_status, ntmp(14)
   real(kind=CFD_solver_real_kind) :: rtmp(3)
   character(len=1) :: bx, by, bz
   character(len=32) :: mode_str, decomp_str, scheme_str
   !-------------------------------------------
   nparameters=0
   rparameters=0.d0

   ntmp = 0
   rtmp = 0.0_CFD_solver_real_kind

   If(my_id .eq. 0) then

      open(30, file='Poisson_solver_UI.in', status='old', action='read', iostat=input_status)
      if (input_status /= 0) go to 900

      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) nx_global,ny_global,nz_global,Dimension,&
                 Boundary(1),Boundary(2),Boundary(3)
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) npx0,npy0,npz0,LAP
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) SLx, SLy, SLz
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) Repeat_count
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) Decomposition
      read(30,*,err=900,end=900)
      read(30,*,err=900,end=900) Scheme_type
      close(30)

   Endif

   if (my_id == 0) then
      ntmp(1)= nx_global
      ntmp(2)= ny_global
      ntmp(3)= nz_global
      ntmp(4)= Dimension
      ntmp(5)= Boundary(1)
      ntmp(6)= Boundary(2)
      ntmp(7)= Boundary(3)
      ntmp(8)= LAP
      ntmp(9)= npx0
      ntmp(10)= npy0
      ntmp(11)= npz0
      ntmp(12)= Decomposition
      ntmp(13)= Scheme_type
      ntmp(14)= Repeat_count


   end if

   call MPI_bcast(ntmp(1),size(ntmp),MPI_INTEGER,0, MPI_COMM_WORLD,ierr)


   nx_global  = ntmp(1)
   ny_global  = ntmp(2)
   nz_global  = ntmp(3)
   Dimension = ntmp(4)
   Boundary(1)= ntmp(5)
   Boundary(2)= ntmp(6)
   Boundary(3)= ntmp(7)
   LAP  = ntmp(8)
   npx0= ntmp(9)
   npy0= ntmp(10)
   npz0= ntmp(11)
   Decomposition    = ntmp(12)
   Scheme_type    = ntmp(13)
   Repeat_count   = ntmp(14)


   if (my_id == 0) then
      rtmp(1)= SLx
      rtmp(2)= SLy
      rtmp(3)= SLz

   end if

   call MPI_bcast(rtmp(1),size(rtmp),CFD_solver_real_mpi_type,0, MPI_COMM_WORLD,ierr)

   SLx= rtmp(1)
   SLy= rtmp(2)
   SLz= rtmp(3)

   call validate_dimension_parameters

   if (my_id == 0) then

      bx = '-'
      by = '-'
      bz = '-'
      if (Boundary(1) == 0) bx = 'P'; if (Boundary(1) == 1) bx = 'D'; if (Boundary(1) == 2) bx = 'N'
      if (Dimension >= 2) then
         if (Boundary(2) == 0) by = 'P'; if (Boundary(2) == 1) by = 'D'; if (Boundary(2) == 2) by = 'N'
      end if
      if (Dimension >= 3) then
         if (Boundary(3) == 0) bz = 'P'; if (Boundary(3) == 1) bz = 'D'; if (Boundary(3) == 2) bz = 'N'
      end if

      if (np_size == 1) then
         mode_str = 'Serial'
         decomp_str = '/'
      else
         mode_str = 'Parallel'
         if (Decomposition == 0) then
            decomp_str = 'Block'
         else if (Decomposition == 1) then
            decomp_str = 'Slab'
         else if (Decomposition == 2) then
            decomp_str = 'Pencil'
         else
            decomp_str = 'Unknown'
            call mpi_finalize(ierr)
            stop
         endif
      endif


      if (Scheme_type == 0) then
         scheme_str = 'Pseudo-Spectral'
      else if (Scheme_type == 1) then
         scheme_str = '2nd Central-Difference'
      else if (Scheme_type == 2) then
         scheme_str = '4th Central-Difference'
      else if (Scheme_type == 3) then
         scheme_str = '4th Compact Central-Difference'
      else if (Scheme_type == 4) then
         scheme_str = '6th Central-Difference'
      else if (Scheme_type == 5) then
         scheme_str = '6th Compact Central-Difference'
      else if (Scheme_type == 6) then
         scheme_str = 'Spectral-like Compact Central-Difference'
      else
         scheme_str = 'Unknown'
         call mpi_finalize(ierr)
         stop
      endif

      print *, "+------------------------------------------------------------------------------+"
      print *, "|             CFD toy solver using ISFFT - Poisson example                     |"
      print *, "+------------------------------------------------------------------------------+"
      print *, "| Parameter             | Value                                                |"
      print *, "+-----------------------+------------------------------------------------------+"
      write(*, '(" | Dimension             | ",I0,T80," |")') Dimension
      write(*, '(" | Grid Size             | (",I0,", ",I0,", ",I0,")",T80," |")') nx_global, ny_global, nz_global
      write(*, '(" | Domain Size           | [0.00, ",F0.2,"] x [0.00, ",F0.2,"] x [0.00, ",F0.2,"]",T80," |")') SLx, SLy, SLz
      write(*, '(" | Boundary Condition    | (",A1,", ",A1,", ",A1,") ", &
         "(P,D,N; - = Inactive)",T80," |")') bx, by, bz
      write(*, '(" | Halos                 | ",I0,T80," |")') LAP
      print *, "+-----------------------+------------------------------------------------------+"
      write(*, '(" | Mode                  | ",A,T80," |")') trim(mode_str)
      write(*, '(" | Total Procs           | ",I0,T80," |")') np_size
      write(*, '(" | Process Topology      | (",I0,", ",I0,", ",I0,")",T80," |")') npx0, npy0, npz0
      print *, "+-----------------------+------------------------------------------------------+"
      write(*, '(" | Decomposition         | ",A,T80," |")') trim(decomp_str)
      write(*, '(" | Scheme                | ",A,T80," |")') trim(scheme_str)
      print *, "+------------------------------------------------------------------------------+"

   endif

   return

900 continue
   write(*,'(A)') 'ERROR: Cannot read Poisson_solver_UI.in. Check the input rows.'
   call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
end subroutine read_parameters

!-----------------------------------------------------------------------------------
subroutine validate_dimension_parameters
   use CFD_solver_parameters
   implicit none

   integer :: ierr, axis
   integer :: grid_size(3)
   logical :: invalid

   invalid = .false.
   grid_size = [nx_global, ny_global, nz_global]

   if (Dimension < 1 .or. Dimension > 3) then
      if (my_id == 0) write(*,'(A,I0,A)') &
         'ERROR: Dimension must be 1, 2, or 3; currently received Dimension = ', Dimension, '.'
      invalid = .true.
   end if

   if (min(nx_global, ny_global, nz_global) < 1) then
      if (my_id == 0) write(*,'(A,3(I0,1X))') &
         'ERROR: All three grid sizes must be positive: ', nx_global, ny_global, nz_global
      invalid = .true.
   end if

   if (min(npx0, npy0, npz0) < 1) then
      if (my_id == 0) write(*,'(A,3(I0,1X))') &
         'ERROR: All three MPI partition counts must be positive: ', npx0, npy0, npz0
      invalid = .true.
   end if

   if (min(SLx, SLy, SLz) <= 0.0_CFD_solver_real_kind) then
      if (my_id == 0) write(*,'(A,3(ES12.4,1X))') &
         'ERROR: All three domain lengths must be positive: ', SLx, SLy, SLz
      invalid = .true.
   end if

   if (Repeat_count < 1) then
      if (my_id == 0) write(*,'(A,I0,A)') &
         'ERROR: Repeat_count must be at least 1. Received ', Repeat_count, '.'
      invalid = .true.
   end if

   if (Scheme_type < 0 .or. Scheme_type > 6) then
      if (my_id == 0) write(*,'(A,I0,A)') &
         'ERROR: Scheme_type must be between 0 and 6; currently received Scheme_type = ', Scheme_type, '.'
      invalid = .true.
   end if

   if (LAP < 0) then
      if (my_id == 0) write(*,'(A,I0,A)') 'ERROR: LAP must be nonnegative; currently received LAP = ', LAP, '.'
      invalid = .true.
   else if (Dimension >= 1 .and. Dimension <= 3) then
      if (LAP >= minval(grid_size(1:Dimension))) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: LAP must be smaller than every active grid size; currently received LAP = ', LAP, '.'
         invalid = .true.
      end if
   end if

   if (Dimension >= 1 .and. Dimension <= 3) then
      do axis = 1, Dimension
         if (Boundary(axis) < 0 .or. Boundary(axis) > 2) then
            if (my_id == 0) write(*,'(A,I0,A,I0,A)') &
               'ERROR: Boundary(', axis, ') must be 0, 1, or 2; currently received ', Boundary(axis), '.'
            invalid = .true.
         else if (Boundary(axis) /= 0 .and. grid_size(axis) < 2) then
            if (my_id == 0) write(*,'(A,I0,A,I0,A)') &
               'ERROR: Active nonperiodic direction ', axis, &
               ' requires at least two grid points; currently received ', grid_size(axis), '.'
            invalid = .true.
         end if
      end do
   end if

   if (Dimension == 1) then
      if (ny_global /= 1 .or. nz_global /= 1) then
         if (my_id == 0) write(*,'(A,2(I0,1X))') &
            'ERROR: Dimension=1 requires ny_global=nz_global=1; currently received ', ny_global, nz_global
         invalid = .true.
      end if
      if (npy0 /= 1 .or. npz0 /= 1) then
         if (my_id == 0) write(*,'(A,2(I0,1X))') &
            'ERROR: Dimension=1 requires npy0=npz0=1; currently received ', npy0, npz0
         invalid = .true.
      end if
      if (np_size /= 1 .or. npx0 /= 1) then
         if (my_id == 0) write(*,'(A,I0,A,I0,A)') &
            'ERROR: Dimension=1 is serial-only and requires np_size=npx0=1; currently received ', &
            np_size, ' and ', npx0, '.'
         invalid = .true.
      end if
      if (Decomposition /= 0) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: Dimension=1 requires Decomposition=0 (serial); currently received ', Decomposition, '.'
         invalid = .true.
      end if
      if (Boundary(1) == 0 .and. mod(nx_global, 2) /= 0) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: The current 1D periodic r2c kernel requires an even nx_global; currently received ', &
            nx_global, '.'
         invalid = .true.
      end if
   else if (Dimension == 2) then
      if (nz_global /= 1) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: Dimension=2 requires nz_global=1; currently received ', nz_global, '.'
         invalid = .true.
      end if
      if (npz0 /= 1) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: Dimension=2 requires npz0=1; currently received ', npz0, '.'
         invalid = .true.
      end if
      if (Decomposition /= 0 .and. Decomposition /= 1) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: Dimension=2 supports only block (0) or slab (1); currently received Decomposition = ', &
            Decomposition, '.'
         invalid = .true.
      end if
      if (npx0 > nx_global .or. npy0 > ny_global) then
         if (my_id == 0) write(*,'(A,2(I0,1X),A,2(I0,1X))') &
            'ERROR: 2D MPI partitions exceed active grid sizes; partitions = ', npx0, npy0, &
            ' grid = ', nx_global, ny_global
         invalid = .true.
      end if
      if (Decomposition == 1 .and. np_size > 1) then
         if (np_size > nx_global .or. np_size > ny_global) then
            if (my_id == 0) write(*,'(A,I0,A,2(I0,1X))') &
               'ERROR: 2D slab requires np_size <= min(nx_global,ny_global); np_size = ', &
               np_size, ' grid = ', nx_global, ny_global
            invalid = .true.
         else if (mod(nx_global, np_size) /= 0 .or. mod(ny_global, np_size) /= 0) then
            if (my_id == 0) write(*,'(A,I0,A,2(I0,1X))') &
               'ERROR: Current 2D slab transpose requires nx_global and ny_global divisible by np_size = ', &
               np_size, '; grid = ', nx_global, ny_global
            invalid = .true.
         end if
      end if
      if (Boundary(1) == 0 .and. mod(nx_global,2) /= 0) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: The current 2D r2c path requires even nx_global when x is periodic; currently received nx_global = ', &
            nx_global, '.'
         invalid = .true.
      else if (Boundary(1) /= 0 .and. Boundary(2) == 0 .and. mod(ny_global,2) /= 0) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: The current 2D r2c path requires even ny_global when y is the first periodic axis; ' // &
            'currently received ny_global = ', &
            ny_global, '.'
         invalid = .true.
      end if
   else if (Dimension == 3) then
      if (Decomposition < 0 .or. Decomposition > 2) then
         if (my_id == 0) write(*,'(A,I0,A)') &
            'ERROR: Dimension=3 requires Decomposition = 0, 1, or 2; currently received Decomposition = ', Decomposition, '.'
         invalid = .true.
      end if
   end if

   if (npx0*npy0*npz0 /= np_size) then
      if (my_id == 0) write(*,'(A,I0,A,I0,A)') &
         'ERROR: npx0*npy0*npz0 = ', npx0*npy0*npz0, &
         ' does not match the MPI process count ', np_size, '.'
      invalid = .true.
   end if

   if (invalid) then
      call mpi_finalize(ierr)
      stop 1
   end if
end subroutine validate_dimension_parameters
