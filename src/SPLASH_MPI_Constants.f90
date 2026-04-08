module SPLASH_MPI_Constants
    use SPLASH_Precision
    implicit none

    integer,save:: my_id,npx,npy,npz,npx0,npy0,npz0, ID_XP1,ID_XM1,ID_YP1,ID_YM1,ID_ZP1,ID_ZM1, &
       MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z, MPI_COMM_XY,MPI_COMM_XZ,MPI_COMM_YZ, SPLASH_Barrier_level, np_size
    integer,save::  i_offset(0:2048),j_offset(0:2048),k_offset(0:2048),i_nn(0:2048),j_nn(0:2048),k_nn(0:2048)

 
 end module SPLASH_MPI_Constants