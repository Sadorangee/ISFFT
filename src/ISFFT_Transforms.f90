!-----------------------------------------------------------------------
! General-purpose ISFFT transforms.
!
! The current research only requires fast Fourier transforms, so this
! module intentionally exposes FFT interfaces only.  DST/DCT interfaces
! may be added later without changing the FFT API defined here.
!
! Contract:
!   * ISFFT_Initialize must be called first.
!   * The initialized global grid size is the transform size.  Arbitrary
!     array lengths are not accepted because the MPI layouts and buffers
!     are built for that grid.
!   * Physical-space arrays contain local owned points only; no halo or
!     ghost cells are accepted.
!   * This module never repartitions user arrays.  For slab/pencil FFTs,
!     callers must first place physical data in the FFT target process
!     layout with ISFFT_Repartitioning_Forward.  After an inverse FFT,
!     callers may restore their original layout explicitly with
!     ISFFT_Repartitioning_Inverse.
!   * Local/global transposes that form part of the FFT algorithm remain
!     internal to this module; they do not restore the caller's original
!     process layout.
!   * Interfaces ending in _inplace overwrite their real/imaginary input
!     arrays and restore the same physical axis order and process layout.
!     They are interface-level in-place transforms and intentionally reuse
!     the module's existing out-of-place transpose buffers internally.
!   * Slab/pencil spectra remain in the native transposed layout.
!       2D slab:          (y, x)
!       3D slab/pencil:   (z, x, y)
!     Serial/block transforms retain the physical axis order.
!   * R2C spectra keep full-size storage, but only kx=0:nx_global/2 is
!     valid and transformed.
!   * Forward FFTs use the existing ISFFT normalization (1/N per
!     transformed axis); inverse FFTs are unscaled.  The public R2C
!     coefficients are adjusted to match the public C2C/standard DFT
!     normalization at every valid half-spectrum mode.
!   * ISFFT_Get_FFT_Local_Sizes returns three-component local shapes.
!     For 1D/2D transforms, unused trailing components are one.
!-----------------------------------------------------------------------
module ISFFT_Transforms
   use ISFFT_Precision
   use ISFFT_Parameters
   use ISFFT_MPI_Constants
   use ISFFT_Buffer
   use ISFFT_FFT_kernel, only: FFT_1D_r2c_kernel => ISFFT_1D_FFT_r2c, &
      FFT_1D_c2c_kernel => ISFFT_1D_FFT_c2c
   use ISFFT_LocalTranspose, only: local_transpose_01_complex, local_transpose_02_complex
   use ISFFT_GlobalTranspose, only: global_transpose_01_complex, global_transpose_10_complex, &
      global_transpose_02_complex, global_transpose_20_complex
   use ISFFT_Repartitioning, only: ISFFT_Get_Target_Layout
   implicit none

   private
   public :: ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_c2c
   public :: ISFFT_2D_FFT_r2c, ISFFT_2D_FFT_c2c
   public :: ISFFT_3D_FFT_r2c, ISFFT_3D_FFT_c2c
   public :: ISFFT_1D_FFT_r2c_inplace, ISFFT_1D_FFT_c2c_inplace
   public :: ISFFT_2D_FFT_r2c_inplace, ISFFT_2D_FFT_c2c_inplace
   public :: ISFFT_3D_FFT_r2c_inplace, ISFFT_3D_FFT_c2c_inplace
   public :: ISFFT_Get_FFT_Local_Sizes

   integer, save :: cached_plan_generation = -1
   integer, save :: physical_local_size(3) = 0
   integer, save :: spectral_local_size(3) = 0
   integer, save :: target_process_grid(3) = 0

contains

   subroutine ISFFT_Get_FFT_Local_Sizes(physical_size, spectral_size)
      integer, intent(out) :: physical_size(3), spectral_size(3)

      call Ensure_Metadata()
      physical_size = physical_local_size
      spectral_size = spectral_local_size
   end subroutine ISFFT_Get_FFT_Local_Sizes

   subroutine ISFFT_1D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:), phi_hat_im(:)
      integer, intent(in) :: sign

      call Validate_Common(1, sign, .true.)
      call Validate_1D_Shapes(phi, phi_hat_re, phi_hat_im)

      if (sign > 0) then
         f_re(:,1,1) = phi
         f_im(:,1,1) = 0.0_ISFFT_REAL_KIND
         call FFT_1D_r2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), +1)
         call Finalize_R2C_Forward(f_re(:,1,1), f_im(:,1,1))
         phi_hat_re = f_re(:,1,1)
         phi_hat_im = f_im(:,1,1)
      else
         f_re(:,1,1) = phi_hat_re
         f_im(:,1,1) = phi_hat_im
         call Prepare_R2C_Inverse(f_re(:,1,1), f_im(:,1,1))
         call FFT_1D_r2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), -1)
         phi = f_re(:,1,1)
      end if
   end subroutine ISFFT_1D_FFT_r2c

   subroutine ISFFT_1D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:), phi_im_inout(:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:), phi_hat_im(:)
      integer, intent(in) :: sign

      call Validate_Common(1, sign, .false.)
      call Validate_1D_C2C_Shapes(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im)

      if (sign > 0) then
         f_re(:,1,1) = phi_re_inout
         f_im(:,1,1) = phi_im_inout
         call FFT_1D_c2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), +1)
         phi_hat_re = f_re(:,1,1)
         phi_hat_im = f_im(:,1,1)
      else
         f_re(:,1,1) = phi_hat_re
         f_im(:,1,1) = phi_hat_im
         call FFT_1D_c2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), -1)
         phi_re_inout = f_re(:,1,1)
         phi_im_inout = f_im(:,1,1)
      end if
   end subroutine ISFFT_1D_FFT_c2c

   subroutine ISFFT_1D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:), phi_im_inout(:)
      integer, intent(in) :: sign

      call Validate_Common(1, sign, .true.)
      call Validate_1D_Inplace_Shapes(phi_re_inout, phi_im_inout)

      if (sign > 0) then
         f_re(:,1,1) = phi_re_inout
         f_im(:,1,1) = 0.0_ISFFT_REAL_KIND
         call FFT_1D_r2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), +1)
         call Finalize_R2C_Forward(f_re(:,1,1), f_im(:,1,1))
         phi_re_inout = f_re(:,1,1)
         phi_im_inout = f_im(:,1,1)
      else
         f_re(:,1,1) = phi_re_inout
         f_im(:,1,1) = phi_im_inout
         call Prepare_R2C_Inverse(f_re(:,1,1), f_im(:,1,1))
         call FFT_1D_r2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), -1)
         phi_re_inout = f_re(:,1,1)
         phi_im_inout = 0.0_ISFFT_REAL_KIND
      end if
   end subroutine ISFFT_1D_FFT_r2c_inplace

   subroutine ISFFT_1D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:), phi_im_inout(:)
      integer, intent(in) :: sign

      call Validate_Common(1, sign, .false.)
      call Validate_1D_Inplace_Shapes(phi_re_inout, phi_im_inout)

      f_re(:,1,1) = phi_re_inout
      f_im(:,1,1) = phi_im_inout
      call FFT_1D_c2c_kernel(nx_global, f_re(:,1,1), f_im(:,1,1), sign)
      phi_re_inout = f_re(:,1,1)
      phi_im_inout = f_im(:,1,1)
   end subroutine ISFFT_1D_FFT_c2c_inplace

   subroutine ISFFT_2D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(:,:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:,:), phi_hat_im(:,:)
      integer, intent(in) :: sign

      call Validate_Common(2, sign, .true.)
      call Validate_2D_Shapes(phi, phi_hat_re, phi_hat_im)
      if (sign > 0) then
         call Load_2D_Real_Physical(phi)
         call Forward_2D(.true.)
         call Store_2D_Spectrum(phi_hat_re, phi_hat_im)
      else
         call Load_2D_Spectrum(phi_hat_re, phi_hat_im)
         call Inverse_2D(.true.)
         call Store_2D_Real_Physical(phi)
      end if
   end subroutine ISFFT_2D_FFT_r2c

   subroutine ISFFT_2D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:), phi_im_inout(:,:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:,:), phi_hat_im(:,:)
      integer, intent(in) :: sign

      call Validate_Common(2, sign, .false.)
      call Validate_2D_C2C_Shapes(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im)
      if (sign > 0) then
         call Load_2D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Forward_2D(.false.)
         call Store_2D_Spectrum(phi_hat_re, phi_hat_im)
      else
         call Load_2D_Spectrum(phi_hat_re, phi_hat_im)
         call Inverse_2D(.false.)
         call Store_2D_Complex_Physical(phi_re_inout, phi_im_inout)
      end if
   end subroutine ISFFT_2D_FFT_c2c

   subroutine ISFFT_2D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:), phi_im_inout(:,:)
      integer, intent(in) :: sign

      call Validate_Common(2, sign, .true.)
      call Validate_2D_Inplace_Shapes(phi_re_inout, phi_im_inout)
      if (sign > 0) then
         call Load_2D_Real_Physical(phi_re_inout)
         call Forward_2D(.true.)
         call Native_To_Physical_Layout_2D()
         call Store_2D_Complex_Physical(phi_re_inout, phi_im_inout)
      else
         call Load_2D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Physical_To_Native_Layout_2D()
         call Inverse_2D(.true.)
         call Store_2D_Real_Physical(phi_re_inout)
         phi_im_inout = 0.0_ISFFT_REAL_KIND
      end if
   end subroutine ISFFT_2D_FFT_r2c_inplace

   subroutine ISFFT_2D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:), phi_im_inout(:,:)
      integer, intent(in) :: sign

      call Validate_Common(2, sign, .false.)
      call Validate_2D_Inplace_Shapes(phi_re_inout, phi_im_inout)
      if (sign > 0) then
         call Load_2D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Forward_2D(.false.)
         call Native_To_Physical_Layout_2D()
      else
         call Load_2D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Physical_To_Native_Layout_2D()
         call Inverse_2D(.false.)
      end if
      call Store_2D_Complex_Physical(phi_re_inout, phi_im_inout)
   end subroutine ISFFT_2D_FFT_c2c_inplace

   subroutine ISFFT_3D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:,:,:), phi_hat_im(:,:,:)
      integer, intent(in) :: sign

      call Validate_Common(3, sign, .true.)
      call Validate_3D_Shapes(phi, phi_hat_re, phi_hat_im)
      if (sign > 0) then
         call Load_3D_Real_Physical(phi)
         call Forward_3D(.true.)
         call Store_3D_Spectrum(phi_hat_re, phi_hat_im)
      else
         call Load_3D_Spectrum(phi_hat_re, phi_hat_im)
         call Inverse_3D(.true.)
         call Store_3D_Real_Physical(phi)
      end if
   end subroutine ISFFT_3D_FFT_r2c

   subroutine ISFFT_3D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_hat_re(:,:,:), phi_hat_im(:,:,:)
      integer, intent(in) :: sign

      call Validate_Common(3, sign, .false.)
      call Validate_3D_C2C_Shapes(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im)
      if (sign > 0) then
         call Load_3D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Forward_3D(.false.)
         call Store_3D_Spectrum(phi_hat_re, phi_hat_im)
      else
         call Load_3D_Spectrum(phi_hat_re, phi_hat_im)
         call Inverse_3D(.false.)
         call Store_3D_Complex_Physical(phi_re_inout, phi_im_inout)
      end if
   end subroutine ISFFT_3D_FFT_c2c

   subroutine ISFFT_3D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      integer, intent(in) :: sign

      call Validate_Common(3, sign, .true.)
      call Validate_3D_Inplace_Shapes(phi_re_inout, phi_im_inout)
      if (sign > 0) then
         call Load_3D_Real_Physical(phi_re_inout)
         call Forward_3D(.true.)
         call Native_To_Physical_Layout_3D()
         call Store_3D_Complex_Physical(phi_re_inout, phi_im_inout)
      else
         call Load_3D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Physical_To_Native_Layout_3D()
         call Inverse_3D(.true.)
         call Store_3D_Real_Physical(phi_re_inout)
         phi_im_inout = 0.0_ISFFT_REAL_KIND
      end if
   end subroutine ISFFT_3D_FFT_r2c_inplace

   subroutine ISFFT_3D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
      integer, intent(in) :: sign

      call Validate_Common(3, sign, .false.)
      call Validate_3D_Inplace_Shapes(phi_re_inout, phi_im_inout)
      if (sign > 0) then
         call Load_3D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Forward_3D(.false.)
         call Native_To_Physical_Layout_3D()
      else
         call Load_3D_Complex_Physical(phi_re_inout, phi_im_inout)
         call Physical_To_Native_Layout_3D()
         call Inverse_3D(.false.)
      end if
      call Store_3D_Complex_Physical(phi_re_inout, phi_im_inout)
   end subroutine ISFFT_3D_FFT_c2c_inplace

   subroutine Ensure_Metadata()
      integer :: target_npx, target_npy, target_npz

      if (ISFFT_plan_generation <= 0 .or. .not. Buffer_3DFFT_Allocated) then
         call Transform_Abort('ISFFT_Initialize must be called before an FFT interface')
      end if
      if (cached_plan_generation == ISFFT_plan_generation) return

      call ISFFT_Get_Target_Layout(target_npx, target_npy, target_npz)
      target_process_grid = [target_npx, target_npy, target_npz]
      physical_local_size = shape(f_re)
      select case (Dimension)
      case (1)
         spectral_local_size = [nx_global, 1, 1]
      case (2)
         if (np_size == 1 .or. If_decomp == 0) then
            spectral_local_size = [size(f_re,1), size(f_re,2), 1]
         else
            spectral_local_size = [size(f_hat_y,1), size(f_hat_y,2), 1]
         end if
      case (3)
         if (np_size == 1 .or. If_decomp == 0) then
            spectral_local_size = shape(f_re)
         else
            spectral_local_size = shape(f_hat)
         end if
      case default
         call Transform_Abort('initialized Dimension must be 1, 2, or 3')
      end select
      cached_plan_generation = ISFFT_plan_generation
   end subroutine Ensure_Metadata

   subroutine Validate_Common(expected_dimension, sign, is_r2c)
      integer, intent(in) :: expected_dimension, sign
      logical, intent(in) :: is_r2c

      call Ensure_Metadata()
      if (Dimension /= expected_dimension) then
         call Transform_Abort('FFT interface rank does not match the initialized Dimension')
      end if
      if (any([npx0,npy0,npz0] /= target_process_grid)) then
         call Transform_Abort('current process layout is not the FFT target layout; call ' // &
            'ISFFT_Repartitioning_Forward on the physical array before this transform')
      end if
      if (sign /= 1 .and. sign /= -1) call Transform_Abort('FFT sign must be +1 or -1')
      if (is_r2c .and. mod(nx_global,2) /= 0) then
         call Transform_Abort('R2C requires an even initialized nx_global')
      end if
      if (expected_dimension == 1 .and. np_size /= 1) then
         call Transform_Abort('the current 1D FFT interface is serial-only')
      end if
   end subroutine Validate_Common

   subroutine Validate_1D_Shapes(phi, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi(:), hat_re(:), hat_im(:)
      if (size(phi) /= physical_local_size(1)) call Transform_Abort('invalid 1D physical array size')
      if (size(hat_re) /= spectral_local_size(1) .or. size(hat_im) /= spectral_local_size(1)) &
         call Transform_Abort('invalid 1D spectral array size')
   end subroutine Validate_1D_Shapes

   subroutine Validate_1D_C2C_Shapes(phi_re_in, phi_im_in, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:), phi_im_in(:), hat_re(:), hat_im(:)
      call Validate_1D_Shapes(phi_re_in, hat_re, hat_im)
      if (size(phi_im_in) /= physical_local_size(1)) call Transform_Abort('invalid 1D imaginary physical array size')
   end subroutine Validate_1D_C2C_Shapes

   subroutine Validate_1D_Inplace_Shapes(phi_re_in, phi_im_in)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:), phi_im_in(:)
      if (size(phi_re_in) /= physical_local_size(1) .or. size(phi_im_in) /= physical_local_size(1)) &
         call Transform_Abort('invalid 1D inplace array size')
   end subroutine Validate_1D_Inplace_Shapes

   subroutine Validate_2D_Shapes(phi, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi(:,:), hat_re(:,:), hat_im(:,:)
      if (any(shape(phi) /= physical_local_size(1:2))) call Transform_Abort('invalid 2D physical array shape')
      if (any(shape(hat_re) /= spectral_local_size(1:2)) .or. &
          any(shape(hat_im) /= spectral_local_size(1:2))) call Transform_Abort('invalid 2D spectral array shape')
   end subroutine Validate_2D_Shapes

   subroutine Validate_2D_C2C_Shapes(phi_re_in, phi_im_in, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:), phi_im_in(:,:), hat_re(:,:), hat_im(:,:)
      call Validate_2D_Shapes(phi_re_in, hat_re, hat_im)
      if (any(shape(phi_im_in) /= physical_local_size(1:2))) &
         call Transform_Abort('invalid 2D imaginary physical array shape')
   end subroutine Validate_2D_C2C_Shapes

   subroutine Validate_2D_Inplace_Shapes(phi_re_in, phi_im_in)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:), phi_im_in(:,:)
      if (any(shape(phi_re_in) /= physical_local_size(1:2)) .or. &
          any(shape(phi_im_in) /= physical_local_size(1:2))) &
         call Transform_Abort('invalid 2D inplace array shape')
   end subroutine Validate_2D_Inplace_Shapes

   subroutine Validate_3D_Shapes(phi, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi(:,:,:), hat_re(:,:,:), hat_im(:,:,:)
      if (any(shape(phi) /= physical_local_size)) call Transform_Abort('invalid 3D physical array shape')
      if (any(shape(hat_re) /= spectral_local_size) .or. any(shape(hat_im) /= spectral_local_size)) &
         call Transform_Abort('invalid 3D spectral array shape')
   end subroutine Validate_3D_Shapes

   subroutine Validate_3D_C2C_Shapes(phi_re_in, phi_im_in, hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
      real(kind=ISFFT_REAL_KIND), intent(in) :: hat_re(:,:,:), hat_im(:,:,:)
      call Validate_3D_Shapes(phi_re_in, hat_re, hat_im)
      if (any(shape(phi_im_in) /= physical_local_size)) &
         call Transform_Abort('invalid 3D imaginary physical array shape')
   end subroutine Validate_3D_C2C_Shapes

   subroutine Validate_3D_Inplace_Shapes(phi_re_in, phi_im_in)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)
      if (any(shape(phi_re_in) /= physical_local_size) .or. any(shape(phi_im_in) /= physical_local_size)) &
         call Transform_Abort('invalid 3D inplace array shape')
   end subroutine Validate_3D_Inplace_Shapes

   subroutine Load_2D_Real_Physical(phi)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi(:,:)

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re(:,:,1) = phi
   end subroutine Load_2D_Real_Physical

   subroutine Load_2D_Complex_Physical(phi_re_in, phi_im_in)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:), phi_im_in(:,:)

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re(:,:,1) = phi_re_in
      f_im(:,:,1) = phi_im_in
   end subroutine Load_2D_Complex_Physical

   subroutine Store_2D_Real_Physical(phi)
      real(kind=ISFFT_REAL_KIND), intent(out) :: phi(:,:)

      phi = f_re(:,:,1)
   end subroutine Store_2D_Real_Physical

   subroutine Store_2D_Complex_Physical(phi_re_out, phi_im_out)
      real(kind=ISFFT_REAL_KIND), intent(out) :: phi_re_out(:,:), phi_im_out(:,:)

      phi_re_out = f_re(:,:,1)
      phi_im_out = f_im(:,:,1)
   end subroutine Store_2D_Complex_Physical

   subroutine Load_3D_Real_Physical(phi)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi(:,:,:)

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re = phi
   end subroutine Load_3D_Real_Physical

   subroutine Load_3D_Complex_Physical(phi_re_in, phi_im_in)
      real(kind=ISFFT_REAL_KIND), intent(in) :: phi_re_in(:,:,:), phi_im_in(:,:,:)

      f_re = 0.0_ISFFT_REAL_KIND
      f_im = 0.0_ISFFT_REAL_KIND
      f_re = phi_re_in
      f_im = phi_im_in
   end subroutine Load_3D_Complex_Physical

   subroutine Store_3D_Real_Physical(phi)
      real(kind=ISFFT_REAL_KIND), intent(out) :: phi(:,:,:)

      phi = f_re
   end subroutine Store_3D_Real_Physical

   subroutine Store_3D_Complex_Physical(phi_re_out, phi_im_out)
      real(kind=ISFFT_REAL_KIND), intent(out) :: phi_re_out(:,:,:), phi_im_out(:,:,:)

      phi_re_out = f_re
      phi_im_out = f_im
   end subroutine Store_3D_Complex_Physical

   subroutine Forward_2D(is_r2c)
      logical, intent(in) :: is_r2c
      if (np_size == 1) then
         call Transform_Serial_Axis(1, +1, is_r2c, f_re, f_im, 2)
         call Transform_Serial_Axis(2, +1, is_r2c, f_re, f_im, 2)
      else if (If_decomp == 0) then
         call Transform_Block_Axis(1, +1, is_r2c, f_re, f_im, 2)
         call Transform_Block_Axis(2, +1, is_r2c, f_re, f_im, 2)
      else
         call Transform_X_Physical(f_re, f_im, +1, is_r2c)
         call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
         call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         call Transform_Y_Layout(f_hat_y, f_hat_im_y, +1, is_r2c)
      end if
   end subroutine Forward_2D

   subroutine Inverse_2D(is_r2c)
      logical, intent(in) :: is_r2c
      if (np_size == 1) then
         call Transform_Serial_Axis(2, -1, is_r2c, f_re, f_im, 2)
         call Transform_Serial_Axis(1, -1, is_r2c, f_re, f_im, 2)
      else if (If_decomp == 0) then
         call Transform_Block_Axis(2, -1, is_r2c, f_re, f_im, 2)
         call Transform_Block_Axis(1, -1, is_r2c, f_re, f_im, 2)
      else
         call Transform_Y_Layout(f_hat_y, f_hat_im_y, -1, is_r2c)
         call global_transpose_10_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
         call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_re, f_im)
         call Transform_X_Physical(f_re, f_im, -1, is_r2c)
      end if
   end subroutine Inverse_2D

   subroutine Physical_To_Native_Layout_2D()
      if (np_size == 1 .or. If_decomp == 0) return

      call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
      call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
   end subroutine Physical_To_Native_Layout_2D

   subroutine Native_To_Physical_Layout_2D()
      if (np_size == 1 .or. If_decomp == 0) return

      call global_transpose_10_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
      call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_re, f_im)
   end subroutine Native_To_Physical_Layout_2D

   subroutine Store_2D_Spectrum(hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(out) :: hat_re(:,:), hat_im(:,:)
      if (np_size == 1 .or. If_decomp == 0) then
         hat_re = f_re(:,:,1)
         hat_im = f_im(:,:,1)
      else
         hat_re = f_hat_y(:,:,1)
         hat_im = f_hat_im_y(:,:,1)
      end if
   end subroutine Store_2D_Spectrum

   subroutine Load_2D_Spectrum(hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: hat_re(:,:), hat_im(:,:)
      if (np_size == 1 .or. If_decomp == 0) then
         f_re(:,:,1) = hat_re
         f_im(:,:,1) = hat_im
      else
         f_hat_y(:,:,1) = hat_re
         f_hat_im_y(:,:,1) = hat_im
      end if
   end subroutine Load_2D_Spectrum

   subroutine Forward_3D(is_r2c)
      logical, intent(in) :: is_r2c
      if (np_size == 1) then
         call Transform_Serial_Axis(1, +1, is_r2c, f_re, f_im, 3)
         call Transform_Serial_Axis(2, +1, is_r2c, f_re, f_im, 3)
         call Transform_Serial_Axis(3, +1, is_r2c, f_re, f_im, 3)
      else if (If_decomp == 0) then
         call Transform_Block_Axis(1, +1, is_r2c, f_re, f_im, 3)
         call Transform_Block_Axis(2, +1, is_r2c, f_re, f_im, 3)
         call Transform_Block_Axis(3, +1, is_r2c, f_re, f_im, 3)
      else
         call Transform_X_Physical(f_re, f_im, +1, is_r2c)
         if (If_decomp == 1) then
            call local_transpose_01_complex(f_re, f_im, f_hat_y, f_hat_im_y)
         else
            call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
            call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         end if
         call Transform_Y_Layout(f_hat_y, f_hat_im_y, +1, is_r2c)
         if (If_decomp == 1) then
            call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
            call global_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat, f_hat_im)
         else
            call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_b, phi_im_transposed_b)
            call global_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, f_hat, f_hat_im)
         end if
         call Transform_Z_Layout(f_hat, f_hat_im, +1, is_r2c)
      end if
   end subroutine Forward_3D

   subroutine Inverse_3D(is_r2c)
      logical, intent(in) :: is_r2c
      if (np_size == 1) then
         call Transform_Serial_Axis(3, -1, is_r2c, f_re, f_im, 3)
         call Transform_Serial_Axis(2, -1, is_r2c, f_re, f_im, 3)
         call Transform_Serial_Axis(1, -1, is_r2c, f_re, f_im, 3)
      else if (If_decomp == 0) then
         call Transform_Block_Axis(3, -1, is_r2c, f_re, f_im, 3)
         call Transform_Block_Axis(2, -1, is_r2c, f_re, f_im, 3)
         call Transform_Block_Axis(1, -1, is_r2c, f_re, f_im, 3)
      else
         call Transform_Z_Layout(f_hat, f_hat_im, -1, is_r2c)
         if (If_decomp == 1) then
            call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_a, phi_im_transposed_a)
            call local_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         else
            call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_b, phi_im_transposed_b)
            call local_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, f_hat_y, f_hat_im_y)
         end if
         call Transform_Y_Layout(f_hat_y, f_hat_im_y, -1, is_r2c)
         if (If_decomp == 1) then
            call local_transpose_01_complex(f_hat_y, f_hat_im_y, f_re, f_im)
         else
            call global_transpose_10_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
            call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_re, f_im)
         end if
         call Transform_X_Physical(f_re, f_im, -1, is_r2c)
      end if
   end subroutine Inverse_3D

   subroutine Physical_To_Native_Layout_3D()
      if (np_size == 1 .or. If_decomp == 0) return

      if (If_decomp == 1) then
         call local_transpose_01_complex(f_re, f_im, f_hat_y, f_hat_im_y)
         call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
         call global_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat, f_hat_im)
      else
         call local_transpose_01_complex(f_re, f_im, phi_re_transposed_a, phi_im_transposed_a)
         call global_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         call local_transpose_02_complex(f_hat_y, f_hat_im_y, phi_re_transposed_b, phi_im_transposed_b)
         call global_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, f_hat, f_hat_im)
      end if
   end subroutine Physical_To_Native_Layout_3D

   subroutine Native_To_Physical_Layout_3D()
      if (np_size == 1 .or. If_decomp == 0) return

      if (If_decomp == 1) then
         call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_a, phi_im_transposed_a)
         call local_transpose_02_complex(phi_re_transposed_a, phi_im_transposed_a, f_hat_y, f_hat_im_y)
         call local_transpose_01_complex(f_hat_y, f_hat_im_y, f_re, f_im)
      else
         call global_transpose_20_complex(f_hat, f_hat_im, phi_re_transposed_b, phi_im_transposed_b)
         call local_transpose_02_complex(phi_re_transposed_b, phi_im_transposed_b, f_hat_y, f_hat_im_y)
         call global_transpose_10_complex(f_hat_y, f_hat_im_y, phi_re_transposed_a, phi_im_transposed_a)
         call local_transpose_01_complex(phi_re_transposed_a, phi_im_transposed_a, f_re, f_im)
      end if
   end subroutine Native_To_Physical_Layout_3D

   subroutine Store_3D_Spectrum(hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(out) :: hat_re(:,:,:), hat_im(:,:,:)
      if (np_size == 1 .or. If_decomp == 0) then
         hat_re = f_re
         hat_im = f_im
      else
         hat_re = f_hat
         hat_im = f_hat_im
      end if
   end subroutine Store_3D_Spectrum

   subroutine Load_3D_Spectrum(hat_re, hat_im)
      real(kind=ISFFT_REAL_KIND), intent(in) :: hat_re(:,:,:), hat_im(:,:,:)
      if (np_size == 1 .or. If_decomp == 0) then
         f_re = hat_re
         f_im = hat_im
      else
         f_hat = hat_re
         f_hat_im = hat_im
      end if
   end subroutine Load_3D_Spectrum

   subroutine Transform_X_Physical(a_re, a_im, sign, is_r2c)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      logical, intent(in) :: is_r2c
      integer :: j, k

      do k = 1, size(a_re,3)
         do j = 1, size(a_re,2)
            if (is_r2c) then
               if (sign < 0) call Prepare_R2C_Inverse(a_re(:,j,k), a_im(:,j,k))
               call FFT_1D_r2c_kernel(nx_global, a_re(:,j,k), a_im(:,j,k), sign)
               if (sign > 0) call Finalize_R2C_Forward(a_re(:,j,k), a_im(:,j,k))
            else
               call FFT_1D_c2c_kernel(nx_global, a_re(:,j,k), a_im(:,j,k), sign)
            end if
         end do
      end do
   end subroutine Transform_X_Physical

   subroutine Transform_Y_Layout(a_re, a_im, sign, is_r2c)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      logical, intent(in) :: is_r2c
      integer :: i, k, i_stop, x_offset

      i_stop = size(a_re,2)
      if (is_r2c) then
         x_offset = npy*size(a_re,2)
         i_stop = max(0, min(size(a_re,2), halfspectrum_x-x_offset))
      end if
      do k = 1, size(a_re,3)
         do i = 1, i_stop
            call FFT_1D_c2c_kernel(ny_global, a_re(:,i,k), a_im(:,i,k), sign)
         end do
      end do
   end subroutine Transform_Y_Layout

   subroutine Transform_Z_Layout(a_re, a_im, sign, is_r2c)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer, intent(in) :: sign
      logical, intent(in) :: is_r2c
      integer :: j, k, j_stop, x_offset

      j_stop = size(a_re,2)
      if (is_r2c) then
         x_offset = npy*size(a_re,2)
         j_stop = max(0, min(size(a_re,2), halfspectrum_x-x_offset))
      end if
      do k = 1, size(a_re,3)
         do j = 1, j_stop
            call FFT_1D_c2c_kernel(nz_global, a_re(:,j,k), a_im(:,j,k), sign)
         end do
      end do
   end subroutine Transform_Z_Layout

   subroutine Transform_Serial_Axis(axis, sign, is_r2c, a_re, a_im, active_dimension)
      integer, intent(in) :: axis, sign, active_dimension
      logical, intent(in) :: is_r2c
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer :: i, j, k, i_stop

      select case (axis)
      case (1)
         do k = 1, merge(nz_global,1,active_dimension == 3)
            do j = 1, merge(ny_global,1,active_dimension >= 2)
               if (is_r2c) then
                  if (sign < 0) call Prepare_R2C_Inverse(a_re(:,j,k), a_im(:,j,k))
                  call FFT_1D_r2c_kernel(nx_global, a_re(:,j,k), a_im(:,j,k), sign)
                  if (sign > 0) call Finalize_R2C_Forward(a_re(:,j,k), a_im(:,j,k))
               else
                  call FFT_1D_c2c_kernel(nx_global, a_re(:,j,k), a_im(:,j,k), sign)
               end if
            end do
         end do
      case (2)
         i_stop = nx_global
         if (is_r2c) i_stop = halfspectrum_x
         do k = 1, merge(nz_global,1,active_dimension == 3)
            do i = 1, i_stop
               call FFT_1D_c2c_kernel(ny_global, a_re(i,:,k), a_im(i,:,k), sign)
            end do
         end do
      case (3)
         i_stop = nx_global
         if (is_r2c) i_stop = halfspectrum_x
         do j = 1, ny_global
            do i = 1, i_stop
               call FFT_1D_c2c_kernel(nz_global, a_re(i,j,:), a_im(i,j,:), sign)
            end do
         end do
      end select
   end subroutine Transform_Serial_Axis

   subroutine Transform_Block_Axis(axis, sign, is_r2c, a_re, a_im, active_dimension)
      integer, intent(in) :: axis, sign, active_dimension
      logical, intent(in) :: is_r2c
      real(kind=ISFFT_REAL_KIND), intent(inout) :: a_re(:,:,:), a_im(:,:,:)
      integer :: i, j, k, gi, gj, gk, ierr, i_stop

      select case (axis)
      case (1)
         do k = 1, merge(nz,1,active_dimension == 3)
            do j = 1, merge(ny,1,active_dimension >= 2)
               x_sendbuf = 0.0_ISFFT_REAL_KIND
               x_im_sendbuf = 0.0_ISFFT_REAL_KIND
               do i = 1, nx
                  gi = i_offset(npx)+i-1
                  x_sendbuf(gi) = a_re(i,j,k)
                  if (.not. is_r2c .or. sign < 0) x_im_sendbuf(gi) = a_im(i,j,k)
               end do
               call MPI_Allreduce(x_sendbuf, x_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
               if (.not. is_r2c .or. sign < 0) then
                  call MPI_Allreduce(x_im_sendbuf, x_im_recvbuf, nx_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_X, ierr)
               else
                  x_im_recvbuf = 0.0_ISFFT_REAL_KIND
               end if
               if (is_r2c) then
                  if (sign < 0) call Prepare_R2C_Inverse(x_recvbuf, x_im_recvbuf)
                  call FFT_1D_r2c_kernel(nx_global, x_recvbuf, x_im_recvbuf, sign)
                  if (sign > 0) call Finalize_R2C_Forward(x_recvbuf, x_im_recvbuf)
               else
                  call FFT_1D_c2c_kernel(nx_global, x_recvbuf, x_im_recvbuf, sign)
               end if
               do i = 1, nx
                  gi = i_offset(npx)+i-1
                  a_re(i,j,k) = x_recvbuf(gi)
                  a_im(i,j,k) = x_im_recvbuf(gi)
               end do
            end do
         end do
      case (2)
         i_stop = nx
         if (is_r2c) i_stop = Active_Local_X_Stop()
         do k = 1, merge(nz,1,active_dimension == 3)
            do i = 1, i_stop
               y_sendbuf = 0.0_ISFFT_REAL_KIND
               y_im_sendbuf = 0.0_ISFFT_REAL_KIND
               do j = 1, ny
                  gj = j_offset(npy)+j-1
                  y_sendbuf(gj) = a_re(i,j,k)
                  y_im_sendbuf(gj) = a_im(i,j,k)
               end do
               call MPI_Allreduce(y_sendbuf, y_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
               call MPI_Allreduce(y_im_sendbuf, y_im_recvbuf, ny_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Y, ierr)
               call FFT_1D_c2c_kernel(ny_global, y_recvbuf, y_im_recvbuf, sign)
               do j = 1, ny
                  gj = j_offset(npy)+j-1
                  a_re(i,j,k) = y_recvbuf(gj)
                  a_im(i,j,k) = y_im_recvbuf(gj)
               end do
            end do
         end do
      case (3)
         i_stop = nx
         if (is_r2c) i_stop = Active_Local_X_Stop()
         do j = 1, ny
            do i = 1, i_stop
               z_sendbuf = 0.0_ISFFT_REAL_KIND
               z_im_sendbuf = 0.0_ISFFT_REAL_KIND
               do k = 1, nz
                  gk = k_offset(npz)+k-1
                  z_sendbuf(gk) = a_re(i,j,k)
                  z_im_sendbuf(gk) = a_im(i,j,k)
               end do
               call MPI_Allreduce(z_sendbuf, z_recvbuf, nz_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Z, ierr)
               call MPI_Allreduce(z_im_sendbuf, z_im_recvbuf, nz_global, ISFFT_DATA_TYPE, MPI_SUM, MPI_COMM_Z, ierr)
               call FFT_1D_c2c_kernel(nz_global, z_recvbuf, z_im_recvbuf, sign)
               do k = 1, nz
                  gk = k_offset(npz)+k-1
                  a_re(i,j,k) = z_recvbuf(gk)
                  a_im(i,j,k) = z_im_recvbuf(gk)
               end do
            end do
         end do
      end select
   end subroutine Transform_Block_Axis

   integer function Active_Local_X_Stop() result(stop_idx)
      stop_idx = max(0, min(nx, halfspectrum_x-i_offset(npx)+1))
   end function Active_Local_X_Stop

   subroutine Finalize_R2C_Forward(line_re, line_im)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(:), line_im(:)
      if (halfspectrum_x > 2) then
         line_re(2:halfspectrum_x-1) = 0.5_ISFFT_REAL_KIND*line_re(2:halfspectrum_x-1)
         line_im(2:halfspectrum_x-1) = 0.5_ISFFT_REAL_KIND*line_im(2:halfspectrum_x-1)
      end if
      if (halfspectrum_x < size(line_re)) then
         line_re(halfspectrum_x+1:) = 0.0_ISFFT_REAL_KIND
         line_im(halfspectrum_x+1:) = 0.0_ISFFT_REAL_KIND
      end if
   end subroutine Finalize_R2C_Forward

   subroutine Prepare_R2C_Inverse(line_re, line_im)
      real(kind=ISFFT_REAL_KIND), intent(inout) :: line_re(:), line_im(:)
      if (halfspectrum_x > 2) then
         line_re(2:halfspectrum_x-1) = 2.0_ISFFT_REAL_KIND*line_re(2:halfspectrum_x-1)
         line_im(2:halfspectrum_x-1) = 2.0_ISFFT_REAL_KIND*line_im(2:halfspectrum_x-1)
      end if
   end subroutine Prepare_R2C_Inverse

   subroutine Transform_Abort(message)
      character(len=*), intent(in) :: message
      integer :: ierr

      if (my_id == 0) write(*,'(A,A)') 'ISFFT transform error: ', trim(message)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      error stop 1
   end subroutine Transform_Abort

end module ISFFT_Transforms
