!!!----------------------------------------
include'Poisson_source.f90'
include'IO_Poisson.f90'
!!!----------------------------------------
subroutine main_solver
   Use CFD_solver_parameters
   use CFD_solver_constants
   use ISFFT
   implicit none

   real(kind=CFD_solver_real_kind), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP):: phi,f
   integer :: Global_Grid_Size(3)
   real(kind=CFD_solver_real_kind) :: Computational_Domain_Size(3)
   integer :: Process_Grid_Size(3)
   integer :: ierr
   real(kind=CFD_solver_real_kind) :: t1, t2, solver_time_local, solver_time_max

   phi = 0.d0
   f = 0.d0

   if (Boundary(1) .eq. 0) hx=SLx/nx_global
   if ((Boundary(1) .eq. 1) .or. (Boundary(1) .eq. 2)) hx=SLx/(nx_global-1)

   hy = 1.0_CFD_solver_real_kind
   if (Dimension >= 2) then
      if (Boundary(2) .eq. 0) hy=SLy/ny_global
      if ((Boundary(2) .eq. 1) .or. (Boundary(2) .eq. 2)) hy=SLy/(ny_global-1)
   end if

   hz = 1.0_CFD_solver_real_kind
   if (Dimension >= 3) then
      if (Boundary(3) .eq. 0) hz=SLz/nz_global
      if ((Boundary(3) .eq. 1) .or. (Boundary(3) .eq. 2)) hz=SLz/(nz_global-1)
   end if

   if (my_id == 0) write(*,'(A,I0)') 'Number of Poisson solver calls: ', Repeat_count

   Global_Grid_Size = [nx_global, ny_global, nz_global]
   Computational_Domain_Size = [SLx, SLy, SLz]
   Process_Grid_Size = [npx0, npy0, npz0]
   
   call ISFFT_Initialize( &
      dimensions=Dimension, global_grid_size=Global_Grid_Size, &
      domain_lengths=Computational_Domain_Size, &
      boundary_types=Boundary, discretization_scheme=Scheme_type, &
      process_grid_size=Process_Grid_Size, decomposition_type=Decomposition, &
      halo_width=LAP, rank=my_id, &
      planning_level=ISFFT_PLAN_MEASURE, &
      prepare_scope=ISFFT_USE_REDISTRIBUTION_POISSON)
   solver_time_local = 0.0_CFD_solver_real_kind

   !!!---------------------------------------------------------------
   !!!-----------------------Main solver-----------------------------
   !!!---------------------------------------------------------------
   do Istep = 1, Repeat_count  ! This toy example uses an explicit call count for clarity.

      ! In real CFD simulations, the source term f is updated at each time step.
      call Poisson_source(f)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      t1 = MPI_Wtime()
      call ISFFT_Poisson_solver(phi, f)
      t2 = MPI_Wtime()
      solver_time_local = solver_time_local+t2-t1

   end do
   !!!---------------------------------------------------------------
   !!!---------------------------------------------------------------
   !!!---------------------------------------------------------------
   call MPI_Reduce(solver_time_local, solver_time_max, 1, CFD_solver_real_mpi_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
   call ISFFT_Finalize()

   if (Dimension == 1) then
      call IO_Poisson_1D(phi)
   else if (Dimension == 2) then
      call IO_Poisson_2D(phi)
   else
      call IO_Poisson_3D(phi)
   end if
   if (my_id .eq. 0) write(*,'(I0,A)') Dimension, 'D Poisson Solve is done!!'

   if (my_id .eq. 0) write(*,*) 'Elapsed time for Poisson solver: ', solver_time_max, ' seconds'


End subroutine main_solver
!!!----------------------------------------
