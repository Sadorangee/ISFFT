module ISFFT_Buffer
   use ISFFT_Precision
   implicit none

   integer, save :: halfspectrum_x = 0
   integer, save :: halfspectrum_y = 0 ! no-need as the order of R2C in ISFFT is x -> y-> z
   integer, save :: halfspectrum_z = 0 ! i.e. halfspectrum only refers to halfspectrum_x in ISFFT

   !------------------------------ Internal Interface Buffers --------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_internal, phi_internal

   !--------------------------------- FFT1D pre-computation ----------------------------------!
   integer, parameter :: FFT_PRECOMPUTE_MAX = 6

   type :: FFT_Precompute_Slot
      integer :: len = 0
      logical :: active = .false.
      integer, allocatable :: rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: twc(:), tws(:)
   end type FFT_Precompute_Slot

   type(FFT_Precompute_Slot), save :: fft_precompute_slots(FFT_PRECOMPUTE_MAX)
   integer, save :: fft_precompute_count = 0
   integer, save :: fft_precompute_next_slot = 1

   !------------------------------------ General FFT Buffer -----------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_re, f_im

   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: sendbuf_3DFFT, sendbuf_im_3DFFT
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: recvbuf_3DFFT, recvbuf_im_3DFFT
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:) :: sendbuf_3DFFT_complex, recvbuf_3DFFT_complex

   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_a, phi_im_transposed_a
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_b, phi_im_transposed_b
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_c, phi_im_transposed_c
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_d, phi_im_transposed_d

   !------------------------------------ Block FFT Buffer -------------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: x_sendbuf, x_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: x_recvbuf, x_im_recvbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: y_sendbuf, y_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: y_recvbuf, y_im_recvbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: z_sendbuf, z_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: z_recvbuf, z_im_recvbuf

   !--------------------------------- Slab & Pencil FFT Buffer --------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat_x, f_hat_im_x
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat_y, f_hat_im_y
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat, f_hat_im
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_hat_y, phi_hat_im_y
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_hat_x, phi_hat_im_x
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re, phi_im

   logical :: Buffer_3DFFT_Allocated = .false.

   !---------------------------- Local Transpose Tuning Parameters ----------------------------!
   integer, save :: tuned_lt12_n1=0, tuned_lt12_n2=0, tuned_lt12_n3=0
   integer, save :: tuned_lt13_n1=0, tuned_lt13_n2=0, tuned_lt13_n3=0
   integer, save :: opt_lt12_bi=32, opt_lt12_bj=32, opt_lt12_bk=4, opt_lt12_variant=0
   integer, save :: opt_lt13_bi=32, opt_lt13_bj=4, opt_lt13_bk=32, opt_lt13_variant=0
   logical, save :: If_lt12_tuned=.false., If_lt13_tuned=.false.

   logical, save :: enable_transpose_tuning = .true.

end module ISFFT_Buffer
