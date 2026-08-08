module ISFFT_Buffer
   use ISFFT_Precision
   implicit none

   integer, save :: halfspectrum_x = 0
   integer, save :: halfspectrum_y = 0 ! no-need as the order of R2C in ISFFT is x -> y-> z
   integer, save :: halfspectrum_z = 0 ! i.e. halfspectrum only refers to halfspectrum_x in ISFFT

   !------------------------------ Internal Interface Buffers --------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_internal, phi_internal

   !--------------------------------- FFT1D pre-computation ----------------------------------!
   integer, parameter :: FFT_MAX_DIMENSIONS = 3
   integer, parameter :: FFT_MAX_PLANS_PER_DIMENSION = 4
   integer, parameter :: FFT_PRECOMPUTE_MAX = FFT_MAX_DIMENSIONS*FFT_MAX_PLANS_PER_DIMENSION

   type :: FFT_Precompute_Slot
      integer :: len = 0
      integer :: last_used = 0
      logical :: active = .false.
      logical :: is_radix2 = .false.
      logical :: uses_fallback = .false.
      logical :: use_bluestein = .false.
      logical :: use_rader = .false.
      logical :: rader_ready = .false.
      logical :: rader_is_radix2 = .false.
      integer :: num_stages = 0
      integer :: bluestein_m = 0
      integer :: rader_n = 0
      integer :: rader_num_stages = 0
      integer, allocatable :: rev(:)
      integer, allocatable :: stage_radix(:), stage_size(:), stage_stride(:)
      integer, allocatable :: rader_perm(:), rader_out_perm(:), rader_stage_radix(:), rader_rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: twc(:), tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: real_twc(:), real_tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: work_r1(:), work_i1(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: work_r2(:), work_i2(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: aux_r(:), aux_i(:)
      integer, allocatable :: bluestein_rev(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_twc(:), bluestein_tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: chirp_c(:), chirp_s(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: kernel_pos_r(:), kernel_pos_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: kernel_neg_r(:), kernel_neg_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_r(:), bluestein_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: bluestein_work_r(:), bluestein_work_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_twc(:), rader_tws(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_kernel_pos_r(:), rader_kernel_pos_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_kernel_neg_r(:), rader_kernel_neg_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_data_r(:), rader_data_i(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_work_r1(:), rader_work_i1(:)
      real(kind=ISFFT_REAL_KIND), allocatable :: rader_work_r2(:), rader_work_i2(:)
   end type FFT_Precompute_Slot

   type(FFT_Precompute_Slot), save :: fft_precompute_slots(FFT_PRECOMPUTE_MAX)
   integer, save :: fft_precompute_count = 0
   integer, save :: fft_precompute_next_slot = 1
   integer, save :: fft_precompute_use_clock = 0

   !------------------------------------ General Buffer -----------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_re, f_im

   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: sendbuf_3DFFT, sendbuf_im_3DFFT
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: recvbuf_3DFFT, recvbuf_im_3DFFT
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:) :: sendbuf_3DFFT_complex, recvbuf_3DFFT_complex

   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_a, phi_im_transposed_a
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re_transposed_b, phi_im_transposed_b

   !------------------------------------ Block Buffer -------------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: x_sendbuf, x_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: x_recvbuf, x_im_recvbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: y_sendbuf, y_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: y_recvbuf, y_im_recvbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: z_sendbuf, z_im_sendbuf
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:) :: z_recvbuf, z_im_recvbuf

   !--------------------------------- Slab & Pencil Buffer --------------------------------!
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat_x, f_hat_im_x
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat_y, f_hat_im_y
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: f_hat, f_hat_im
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_hat_y, phi_hat_im_y
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_hat_x, phi_hat_im_x
   real(kind=ISFFT_REAL_KIND), allocatable, dimension(:,:,:) :: phi_re, phi_im
   
   !-------------------------- Local Transpose Tuning Parameters --------------------------!
   integer, save :: tuned_lt01_n1=0, tuned_lt01_n2=0, tuned_lt01_n3=0
   integer, save :: tuned_lt02_n1=0, tuned_lt02_n2=0, tuned_lt02_n3=0
   integer, save :: opt_lt01_bi=32, opt_lt01_bj=32, opt_lt01_bk=4, opt_lt01_variant=0
   integer, save :: opt_lt02_bi=32, opt_lt02_bj=4, opt_lt02_bk=32, opt_lt02_variant=0
   logical, save :: If_lt01_tuned=.false., If_lt02_tuned=.false.

   logical, save :: enable_transpose_tuning = .true.

   logical :: Buffer_3DFFT_Allocated = .false.
   logical :: Buffer_3DDST1_Allocated = .false.
   logical :: Buffer_3DDCT1_Allocated = .false.

   logical, save :: enable_FFT = .false.
   logical, save :: enable_DST1 = .false.
   logical, save :: enable_DCT1 = .false.

end module ISFFT_Buffer
