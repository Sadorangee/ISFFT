!------------------------------------------------------------------------------
! ISFFT_Precision
! Fortran integer and real kinds with matching MPI element types.
!------------------------------------------------------------------------------
module ISFFT_Precision
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_INTEGER
   implicit none

   private
   public :: ISFFT_REAL_KIND, ISFFT_INT64_KIND
   public :: ISFFT_REAL_MPI_TYPE, ISFFT_COMPLEX_MPI_TYPE, ISFFT_INTEGER_MPI_TYPE

   integer, parameter :: ISFFT_REAL_KIND = kind(1.0d0)
   !! Wide counts and displacements use the matching MPI_INTEGER8 communication type.
   integer, parameter :: ISFFT_INT64_KIND = selected_int_kind(18)

   integer, parameter :: ISFFT_REAL_MPI_TYPE = MPI_DOUBLE_PRECISION
   integer, parameter :: ISFFT_COMPLEX_MPI_TYPE = MPI_DOUBLE_COMPLEX
   integer, parameter :: ISFFT_INTEGER_MPI_TYPE = MPI_INTEGER

end module ISFFT_Precision
