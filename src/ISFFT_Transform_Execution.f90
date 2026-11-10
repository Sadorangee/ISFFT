!------------------------------------------------------------------------------
! ISFFT_Transform_Execution
! Complete forward/inverse execution and public FFT entry points.
!------------------------------------------------------------------------------
module ISFFT_Transform_Execution
   use mpi
   use ISFFT_Communicator, only: ISFFT_Comm
   use ISFFT_Precision, only: ISFFT_INT64_KIND, ISFFT_REAL_KIND
   use ISFFT_Plans, only: Axis_Plan, Kernel_Plan, public_export, public_import, &
      public_plans, Transform_Plan
   use ISFFT_Problem_State, only: dimensions, internal_layout, Layout, my_id, &
      npx0, npy0, npz0, nx_global, &
      public_fft_modes, public_layout, topology_ready
   use ISFFT_Workspace, only: Get_Complex, Load_Field, Park_Field, poisson_field, &
      Spectral_Field
   use ISFFT_MPI_Part, only: Axis_Position, Layout_Box, Layout_Shape, Same_Layout
   use ISFFT_Axis_Execution, only: Axis_Execute, Execute_Kernel, Execute_Kernel_Forward
   use ISFFT_Communication, only: Begin_Block_Transfer, Block_Transfer, End_Block_Transfer, Gather_Block_Tile, &
      Move_State, Store_Block_Tile
   implicit none
   private
   public :: Block_Sample_Lines, Can_Direct_Input, Can_Direct_Output, Execute_Block_Axis, Execute_Forward
   public :: Execute_Inverse, ISFFT_1D_FFT_c2c, ISFFT_1D_FFT_c2c_inplace, ISFFT_1D_FFT_r2c, ISFFT_1D_FFT_r2c_inplace
   public :: ISFFT_2D_FFT_c2c, ISFFT_2D_FFT_c2c_inplace, ISFFT_2D_FFT_r2c, ISFFT_2D_FFT_r2c_inplace, ISFFT_3D_FFT_c2c
   public :: ISFFT_3D_FFT_c2c_inplace, ISFFT_3D_FFT_r2c, ISFFT_3D_FFT_r2c_inplace, ISFFT_Get_FFT_Layout, ISFFT_Get_FFT_Local_Sizes

contains

   !------------------------------------------------------------------------------
   !! Public physical and spectral layout queries
   !------------------------------------------------------------------------------
   subroutine ISFFT_Get_FFT_Local_Sizes(physical_size, spectral_size, is_r2c)
      integer, intent(out) :: physical_size(3), spectral_size(3)
      logical, intent(in), optional :: is_r2c
      integer :: which

      if (.not. topology_ready) error stop 'ISFFT must be initialized'
      which = 1
      if (present(is_r2c)) then
         if (is_r2c) which = 2
      end if
      if (.not. public_plans(which)%ready) error stop 'Requested FFT usage was not included in ISFFT_Initialize'
      physical_size = Layout_Shape(internal_layout, my_id)
      if (.not. btest(public_fft_modes(which), 0)) error stop 'Internal FFT output mode was not prepared'
      spectral_size = Layout_Shape(public_layout(which, 1), my_id)
   end subroutine ! ISFFT_Get_FFT_Local_Sizes

   !! global_offset is zero-based. Sizes and offsets use memory-axis order.
   !! axis_order maps memory axes to physical axes.
   subroutine ISFFT_Get_FFT_Layout( &
      local_size, global_offset, &
      axis_order, process_grid, half_axis, &
      is_r2c, original)
      integer, intent(out) :: local_size(3), global_offset(3), axis_order(3), process_grid(3), half_axis
      logical, intent(in) :: is_r2c
      logical, intent(in), optional :: original
      type(Layout) :: l
      integer :: which, mode, n(3), off(3)

      which = 1
      if (is_r2c) which = 2
      mode = 1
      if (present(original)) then
         if (original) mode = 2
      end if
      if (.not. public_plans(which)%ready .or. .not. btest(public_fft_modes(which), mode - 1)) &
         error stop 'FFT mode not prepared'
      l = public_layout(which, mode)
      call Layout_Box(l, my_id, n, off)
      local_size = n(l%axis_order)
      global_offset = off(l%axis_order)
      axis_order = l%axis_order
      process_grid = l%process_grid
      half_axis = 0
      if (is_r2c) half_axis = public_plans(which)%half_axis
   end subroutine ! ISFFT_Get_FFT_Layout

   !------------------------------------------------------------------------------
   !! Public coefficient convention and execution dispatch
   !------------------------------------------------------------------------------
   !! Import and export factors are reciprocal. Keep DC and Nyquist unpaired.
   !! The historical real kernel stores paired interior modes with weight two.
   !! Public FFT coefficients use the standard normalized DFT convention.
   !! Fuse the public/internal R2C convention conversion with the required
   !! split-real interface copy. The shared FFT executor keeps one convention.
   subroutine Copy_Public_R2C(z, re, im, field_layout, exporting)
      complex(ISFFT_REAL_KIND), intent(inout) :: z(:, :, :)
      real(ISFFT_REAL_KIND), intent(inout) :: re(:, :, :), im(:, :, :)
      type(Layout), intent(in) :: field_layout
      logical, intent(in) :: exporting
      integer :: n(3), off(3), i, j, k, axis, mode
      real(ISFFT_REAL_KIND) :: factor, paired_factor

      call Layout_Box(field_layout, my_id, n, off)
      axis = Axis_Position(field_layout%axis_order, 1)
      paired_factor = 2.0_ISFFT_REAL_KIND
      if (exporting) paired_factor = 0.5_ISFFT_REAL_KIND
      do k = 1, size(z, 3)
         do j = 1, size(z, 2)
            do i = 1, size(z, 1)
               select case (axis)
               case (1)
                  mode = off(1) + i - 1
               case (2)
                  mode = off(1) + j - 1
               case (3)
                  mode = off(1) + k - 1
               end select
               factor = 1.0_ISFFT_REAL_KIND
               if (mode > 0 .and. 2 * mode < nx_global) factor = paired_factor
               if (exporting) then
                  re(i, j, k) = factor * real(z(i, j, k), ISFFT_REAL_KIND)
                  im(i, j, k) = factor * aimag(z(i, j, k))
               else
                  z(i, j, k) = factor * cmplx(re(i, j, k), im(i, j, k), ISFFT_REAL_KIND)
               end if
            end do
         end do
      end do
   end subroutine ! Copy_Public_R2C

   !! sign=+1 selects forward and sign=-1 selects inverse.
   !! inplace selects original-ownership output. Query the selected spectral layout.
   subroutine Public_Execute(re, im, hr, hi, sign, is_r2c, inplace)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous :: re(:, :, :), im(:, :, :), hr(:, :, :), hi(:, :, :)
      integer, intent(in) :: sign
      logical, intent(in) :: is_r2c, inplace
      logical :: direct_input, direct_output
      type(Spectral_Field), pointer :: f
      integer :: which, mode, ns(3), phys(3), full_shape(3)
      type(Layout) :: full_layout

      f => poisson_field
      which = 1
      if (is_r2c) which = 2
      if (.not. public_plans(which)%ready) error stop 'ISFFT: FFT plan missing; call ISFFT_Initialize first'
      if (sign /= 1 .and. sign /= -1) error stop 'FFT sign must be +1 or -1'
      if (any([npx0, npy0, npz0] /= internal_layout%process_grid)) error stop 'FFT requires target process layout'
      mode = 1
      if (inplace) mode = 2
      if (.not. btest(public_fft_modes(which), mode - 1)) error stop 'FFT output mode not prepared'
      phys = Layout_Shape(internal_layout, my_id)
      if (any(shape(re) /= phys)) error stop 'FFT physical shape mismatch'
      if (.not. is_r2c .or. inplace) then
         if (any(shape(im) /= phys)) error stop 'FFT imaginary shape mismatch'
      end if
      ns = Layout_Shape(public_layout(which, mode), my_id)
      if (.not. inplace) then
         full_layout = public_layout(which, mode)
         full_layout%global_shape = internal_layout%global_shape
         full_shape = Layout_Shape(full_layout, my_id)
         if ((any(shape(hr) /= ns) .and. any(shape(hr) /= full_shape)) .or. any(shape(hi) /= shape(hr))) &
            error stop 'FFT spectrum must have the queried half or full shape'
      end if

      !! Forward execution and selected public spectrum export
      if (sign > 0) then
         !! Out-of-place forward transforms preserve the caller's physical input.
         direct_input = is_r2c .and. Can_Direct_Input(public_plans(which))
         if (direct_input) then
            call Execute_Forward(public_plans(which), f, re)
         else
            if (is_r2c) then
               call Load_Field(f, re)
            else
               call Get_Complex(f, phys, f%z)
               f%z = cmplx(re, im, ISFFT_REAL_KIND)
            end if
            call Execute_Forward(public_plans(which), f)
         end if
         if (mode == 2) call Move_State(public_export(which, mode), f)
         if (inplace) then
            re = 0
            im = 0
            if (is_r2c) then
               call Copy_Public_R2C(f%z, re, im, public_layout(which, mode), .true.)
            else
               re(1:ns(1), 1:ns(2), 1:ns(3)) = real(f%z, ISFFT_REAL_KIND)
               im(1:ns(1), 1:ns(2), 1:ns(3)) = aimag(f%z)
            end if
         else
            hr = 0
            hi = 0
            if (is_r2c) then
               call Copy_Public_R2C(f%z, hr, hi, public_layout(which, mode), .true.)
            else
               hr(1:ns(1), 1:ns(2), 1:ns(3)) = real(f%z, ISFFT_REAL_KIND)
               hi(1:ns(1), 1:ns(2), 1:ns(3)) = aimag(f%z)
            end if
         end if
      else

         !! Inverse spectrum import and physical output restoration
         call Get_Complex(f, ns, f%z)
         if (inplace) then
            if (is_r2c) then
               call Copy_Public_R2C(f%z, re, im, public_layout(which, mode), .false.)
            else
               f%z = cmplx(re(1:ns(1), 1:ns(2), 1:ns(3)), im(1:ns(1), 1:ns(2), 1:ns(3)), ISFFT_REAL_KIND)
            end if
         else
            if (is_r2c) then
               call Copy_Public_R2C(f%z, hr, hi, public_layout(which, mode), .false.)
            else
               f%z = cmplx(hr(1:ns(1), 1:ns(2), 1:ns(3)), hi(1:ns(1), 1:ns(2), 1:ns(3)), ISFFT_REAL_KIND)
            end if
         end if
         if (mode == 2) call Move_State(public_import(which, mode), f)
         direct_output = is_r2c .and. Can_Direct_Output(public_plans(which))
         if (direct_output) then
            call Execute_Inverse(public_plans(which), f, re)
         else
            call Execute_Inverse(public_plans(which), f)
         end if
         if (is_r2c) then
            if (.not. direct_output) re = f%r
            if (inplace) im = 0
         else
            re = real(f%z, ISFFT_REAL_KIND)
            im = aimag(f%z)
         end if
      end if

      !! Return execution storage to the prepared field pool
      call Park_Field(f)
   end subroutine ! Public_Execute

   !------------------------------------------------------------------------------
   !! One-dimensional public FFT adapters
   !------------------------------------------------------------------------------
   subroutine ISFFT_1D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi(:), phi_hat_re(:), phi_hat_im(:)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 1) error stop 'FFT dimension mismatch'
      p0(1:size(phi, 1), 1:1, 1:1) => phi
      p1(1:size(phi_hat_re, 1), 1:1, 1:1) => phi_hat_re
      p2(1:size(phi_hat_im, 1), 1:1, 1:1) => phi_hat_im
      call Public_Execute(p0, dummy_im, p1, p2, sign, .true., .false.)
   end subroutine ! ISFFT_1D_FFT_r2c

   subroutine ISFFT_1D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:), phi_im_inout(:)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 1) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:1, 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:1, 1:1) => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .true., .true.)
   end subroutine ! ISFFT_1D_FFT_r2c_inplace

   subroutine ISFFT_1D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:), phi_im_inout(:), &
         phi_hat_re(:), &
         phi_hat_im(:)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p3(:, :, :)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 1) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:1, 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:1, 1:1) => phi_im_inout
      p2(1:size(phi_hat_re, 1), 1:1, 1:1) => phi_hat_re
      p3(1:size(phi_hat_im, 1), 1:1, 1:1) => phi_hat_im
      call Public_Execute(p0, p1, p2, p3, sign, .false., .false.)
   end subroutine ! ISFFT_1D_FFT_c2c

   subroutine ISFFT_1D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:), phi_im_inout(:)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 1) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:1, 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:1, 1:1) => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .false., .true.)
   end subroutine ! ISFFT_1D_FFT_c2c_inplace

   !------------------------------------------------------------------------------
   !! Two-dimensional public FFT adapters
   !------------------------------------------------------------------------------
   subroutine ISFFT_2D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi(:, :), phi_hat_re(:, :), phi_hat_im(:, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 2) error stop 'FFT dimension mismatch'
      p0(1:size(phi, 1), 1:size(phi, 2), 1:1) => phi
      p1(1:size(phi_hat_re, 1), 1:size(phi_hat_re, 2), 1:1) => phi_hat_re
      p2(1:size(phi_hat_im, 1), 1:size(phi_hat_im, 2), 1:1) => phi_hat_im
      call Public_Execute(p0, dummy_im, p1, p2, sign, .true., .false.)
   end subroutine ! ISFFT_2D_FFT_r2c

   subroutine ISFFT_2D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :), phi_im_inout(:, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 2) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:size(phi_re_inout, 2), 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:size(phi_im_inout, 2), 1:1) => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .true., .true.)
   end subroutine ! ISFFT_2D_FFT_r2c_inplace

   subroutine ISFFT_2D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :), phi_im_inout(:, :), &
         phi_hat_re(:, :), phi_hat_im(:, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p3(:, :, :)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 2) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:size(phi_re_inout, 2), 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:size(phi_im_inout, 2), 1:1) => phi_im_inout
      p2(1:size(phi_hat_re, 1), 1:size(phi_hat_re, 2), 1:1) => phi_hat_re
      p3(1:size(phi_hat_im, 1), 1:size(phi_hat_im, 2), 1:1) => phi_hat_im
      call Public_Execute(p0, p1, p2, p3, sign, .false., .false.)
   end subroutine ! ISFFT_2D_FFT_c2c

   subroutine ISFFT_2D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :), phi_im_inout(:, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 2) error stop 'FFT dimension mismatch'
      p0(1:size(phi_re_inout, 1), 1:size(phi_re_inout, 2), 1:1) => phi_re_inout
      p1(1:size(phi_im_inout, 1), 1:size(phi_im_inout, 2), 1:1) => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .false., .true.)
   end subroutine ! ISFFT_2D_FFT_c2c_inplace

   !------------------------------------------------------------------------------
   !! Three-dimensional public FFT adapters
   !------------------------------------------------------------------------------
   subroutine ISFFT_3D_FFT_r2c(phi, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi(:, :, :), phi_hat_re(:, :, :), &
         phi_hat_im(:, :, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 3) error stop 'FFT dimension mismatch'
      p0 => phi
      p1 => phi_hat_re
      p2 => phi_hat_im
      call Public_Execute(p0, dummy_im, p1, p2, sign, .true., .false.)
   end subroutine ! ISFFT_3D_FFT_r2c

   subroutine ISFFT_3D_FFT_r2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :, :), phi_im_inout(:, :, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 3) error stop 'FFT dimension mismatch'
      p0 => phi_re_inout
      p1 => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .true., .true.)
   end subroutine ! ISFFT_3D_FFT_r2c_inplace

   subroutine ISFFT_3D_FFT_c2c(phi_re_inout, phi_im_inout, phi_hat_re, phi_hat_im, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :, :), phi_im_inout(:, :, :), &
         phi_hat_re(:, :, :), phi_hat_im(:, :, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p2(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p3(:, :, :)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 3) error stop 'FFT dimension mismatch'
      p0 => phi_re_inout
      p1 => phi_im_inout
      p2 => phi_hat_re
      p3 => phi_hat_im
      call Public_Execute(p0, p1, p2, p3, sign, .false., .false.)
   end subroutine ! ISFFT_3D_FFT_c2c

   subroutine ISFFT_3D_FFT_c2c_inplace(phi_re_inout, phi_im_inout, sign)
      real(kind=ISFFT_REAL_KIND), intent(inout), contiguous, target :: phi_re_inout(:, :, :), phi_im_inout(:, :, :)
      integer, intent(in) :: sign
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p0(:, :, :)
      real(kind=ISFFT_REAL_KIND), pointer, contiguous :: p1(:, :, :)
      real(kind=ISFFT_REAL_KIND) :: dummy_im(0, 0, 0), dummy_re(0, 0, 0)

      if (.not. topology_ready) error stop 'ISFFT: call ISFFT_Initialize before FFT execution'
      if (dimensions /= 3) error stop 'FFT dimension mismatch'
      p0 => phi_re_inout
      p1 => phi_im_inout
      call Public_Execute(p0, p1, dummy_re, dummy_im, sign, .false., .true.)
   end subroutine ! ISFFT_3D_FFT_c2c_inplace

   !------------------------------------------------------------------------------
   !! Distributed block tiles and execution samples
   !------------------------------------------------------------------------------
   !! Sampling visits complete execution tiles, including the actual final tile.
   !! It never changes the physical array strides or invents a new tail width.
   pure integer function Block_Sample_Lines(line_count, batch_width, sample_line_budget) result(count)
      integer, intent(in) :: line_count, batch_width, sample_line_budget
      integer :: tiles, visits, cap

      cap = max(1, batch_width)
      tiles = 0
      if (line_count > 0) tiles = 1 + (line_count - 1) / cap
      visits = tiles
      if (sample_line_budget > 0) visits = min(tiles, max(2, (sample_line_budget + cap - 1) / cap))
      count = line_count
      if (visits < tiles) count = (visits - 1) * cap + line_count - (tiles - 1) * cap
   end function ! Block_Sample_Lines

   !! All gather ranks must visit the same tile schedule, including the final tile.
   !! Tail width and its prepared kernel must agree with execution.
   subroutine Execute_Block_Axis(plan, field, sign, sample_lines)
      type(Axis_Plan), intent(in) :: plan
      type(Spectral_Field), intent(inout) :: field
      integer, intent(in) :: sign
      integer, intent(in), optional :: sample_lines
      type(Block_Transfer) :: transfer
      type(Kernel_Plan) :: kernel
      integer :: lines, capacity, tiles, visits, batch, first, width

      call Begin_Block_Transfer(plan, field, sign, transfer)
      lines = transfer%lines
      capacity = transfer%capacity
      tiles = 0
      if (lines > 0) tiles = 1 + (lines - 1) / capacity
      visits = tiles
      if (present(sample_lines)) then
         if (sample_lines > 0) visits = min(tiles, max(2, (sample_lines + capacity - 1) / capacity))
      end if

      !! Tile loop uses the same geometry for sampling and execution
      do batch = 0, visits - 1
         first = batch * capacity
         if (visits < tiles) first = int(int(batch, ISFFT_INT64_KIND) * (tiles - 1) / (visits - 1)) * capacity
         width = min(capacity, lines - first)
         call Gather_Block_Tile(transfer, field, first, width)
         associate (source => transfer%source, tile_r => transfer%tile_r, tile_z => transfer%tile_z)
            !! Execute the prepared full-width or tail kernel
            kernel = plan%kernel
            if (width /= capacity .and. plan%tail_kernel%transform_length > 0) kernel = plan%tail_kernel
            if (plan%transform_kind == 1) then
               call Execute_Kernel(kernel, tile_r(:, 1:width, :), tile_z(:, 1:width, :), sign)
            else if (source%complex_data) then
               call Execute_Kernel(kernel, tile_z(:, 1:width, :), sign)
            else
               call Execute_Kernel(kernel, tile_r(:, 1:width, :), sign)
            end if
   
         end associate
         call Store_Block_Tile(transfer, first, width)
      end do
      call End_Block_Transfer(transfer, field)
   end subroutine ! Execute_Block_Axis

   !------------------------------------------------------------------------------
   !! Forward and inverse transform pipelines
   !------------------------------------------------------------------------------
   logical function Can_Direct_Input(plan) result(ok)
      type(Transform_Plan), intent(in) :: plan

      ok = plan%axes(1)%transform_kind == 1 .and. .not. plan%axes(1)%tiled_block
      ok = ok .and. Same_Layout(plan%axes(1)%before%source, plan%axes(1)%before%destination)
      ok = ok .and. Same_Layout(plan%axes(1)%pack%source, plan%axes(1)%pack%destination)
   end function ! Can_Direct_Input

   subroutine Execute_Forward(plan, field, input)
      type(Transform_Plan), intent(inout) :: plan
      type(Spectral_Field), intent(inout) :: field
      real(kind=ISFFT_REAL_KIND), optional, contiguous, target, intent(in) :: input(:, :, :)
      complex(kind=ISFFT_REAL_KIND), allocatable :: spectrum(:, :, :)
      integer :: a, n(3), rank, ierr

      if (.not. plan%ready) error stop 'transform plan is not ready'
      if (present(input)) then
         if (.not. Can_Direct_Input(plan)) error stop 'plan cannot consume real input directly'
         call MPI_Comm_rank(ISFFT_Comm, rank, ierr)
         if (any(shape(input) /= Layout_Shape(plan%axes(1)%input, rank))) error stop 'direct input shape'
         call Park_Field(field)
      end if
      do a = 1, plan%dimensions
         if (a == 1 .and. present(input)) then
            n = Layout_Shape(plan%axes(1)%output, rank)
            call Get_Complex(field, n, spectrum)
            call Execute_Kernel_Forward(plan%axes(1)%kernel, input, spectrum)
            call move_alloc(spectrum, field%z)
            call Move_State(plan%axes(1)%unpack, field)
            call Move_State(plan%axes(1)%after, field)
            cycle
         end if
         if (plan%axes(a)%tiled_block) then
            call Execute_Block_Axis(plan%axes(a), field, 1)
            cycle
         end if
         call Move_State(plan%axes(a)%before, field)
         call Move_State(plan%axes(a)%pack, field)
         call Axis_Execute(plan%axes(a), field, 1)
         call Move_State(plan%axes(a)%unpack, field)
         call Move_State(plan%axes(a)%after, field)
      end do
   end subroutine ! Execute_Forward

   logical function Can_Direct_Output(plan) result(ok)
      type(Transform_Plan), intent(in) :: plan

      ok = plan%axes(1)%transform_kind == 1 .and. .not. plan%axes(1)%tiled_block
      ok = ok .and. Same_Layout(plan%axes(1)%reverse_unpack%source, plan%axes(1)%reverse_unpack%destination)
      ok = ok .and. Same_Layout(plan%axes(1)%reverse_after%source, plan%axes(1)%reverse_after%destination)
   end function ! Can_Direct_Output

   !! Reverse both axis order and layout transitions from Execute_Forward.
   subroutine Execute_Inverse(plan, field, output)
      type(Transform_Plan), intent(inout) :: plan
      type(Spectral_Field), intent(inout) :: field
      real(kind=ISFFT_REAL_KIND), optional, contiguous, intent(inout) :: output(:, :, :)
      integer :: a

      if (.not. plan%ready) error stop 'transform plan is not ready'
      do a = plan%dimensions, 1, -1
         if (plan%axes(a)%tiled_block) then
            call Execute_Block_Axis(plan%axes(a), field, -1)
            cycle
         end if
         call Move_State(plan%axes(a)%reverse_before, field)
         call Move_State(plan%axes(a)%reverse_pack, field)
         if (a == 1 .and. present(output)) then
            if (.not. Can_Direct_Output(plan)) error stop 'plan cannot write directly to real output'
            call Execute_Kernel(plan%axes(a)%kernel, output, field%z, -1)
            call Park_Field(field)
            return
         end if
         call Axis_Execute(plan%axes(a), field, -1)
         call Move_State(plan%axes(a)%reverse_unpack, field)
         call Move_State(plan%axes(a)%reverse_after, field)
      end do
   end subroutine ! Execute_Inverse

end module ISFFT_Transform_Execution
