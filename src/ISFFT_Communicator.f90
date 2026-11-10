!------------------------------------------------------------------------------
! ISFFT_Communicator
! Library communicator selection, ownership and release.
!------------------------------------------------------------------------------
! The library owns a duplicate only when a caller supplies comm.
module ISFFT_Communicator
   use mpi
   implicit none
   private
   public :: ISFFT_Comm_Process, ISFFT_Comm, ISFFT_Set_Communicator, ISFFT_Release_Communicator
   integer, save :: ISFFT_Comm = MPI_COMM_WORLD
   integer, save :: ISFFT_Comm_Process = -1
   logical, save :: owned = .false.
contains

   !! Call after MPI_Init. A supplied communicator is duplicated collectively.
   !! Comm_free applies only to the owned duplicate.
   subroutine ISFFT_Set_Communicator(comm)
      integer, intent(in), optional :: comm
      integer :: ierr
      logical :: inter

      if (owned) error stop 'ISFFT communicator already initialized'
      ISFFT_Comm = MPI_COMM_WORLD
      if (present(comm)) then
         !! World rank is used only as a unique diagnostic file identifier.
         call MPI_Comm_rank(MPI_COMM_WORLD, ISFFT_Comm_Process, ierr)
         if (comm == MPI_COMM_NULL) error stop 'ISFFT requires a valid intracommunicator'
         call MPI_Comm_test_inter(comm, inter, ierr)
         if (inter) error stop 'ISFFT requires an intracommunicator'
         call MPI_Comm_dup(comm, ISFFT_Comm, ierr)
         if (ierr /= MPI_SUCCESS) error stop 'ISFFT communicator duplication failed'
         owned = .true.
      end if
   end subroutine ! ISFFT_Set_Communicator

   subroutine ISFFT_Release_Communicator()
      integer :: ierr

      if (owned) call MPI_Comm_free(ISFFT_Comm, ierr)
      owned = .false.
      ISFFT_Comm = MPI_COMM_WORLD
      ISFFT_Comm_Process = -1
   end subroutine ! ISFFT_Release_Communicator
end module ! ISFFT_Communicator
