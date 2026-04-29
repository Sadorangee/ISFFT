module ISFFT_Poisson_Algebra
    use ISFFT_Parameters
    use ISFFT_MPI_Constants
    use ISFFT_Buffer
    implicit none
 
 contains

 !----------------------------------------------------------------------
 subroutine Algebra_Periodic_Poisson_r2c(phi_re_inout, phi_im_inout)

    integer :: i, j, k, ierr
    integer :: kx, ky, kz, k1, k2, k3, k_xi, k_eta, k_zeta, nx0, ny0, nz0

    real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
    real(kind=ISFFT_REAL_KIND) :: sxx, syy, szz
    real(kind=ISFFT_REAL_KIND) :: kxi_scaled, keta_scaled, kzeta_scaled, k_squared_total, inv_k_squared
    real(kind=ISFFT_REAL_KIND), parameter :: PI = 3.1415926535897932384626433832795_ISFFT_REAL_KIND

    ! Pre-compute constants
    sxx = 2.0_ISFFT_REAL_KIND * PI / SLx
    syy = 2.0_ISFFT_REAL_KIND * PI / SLy
    szz = 2.0_ISFFT_REAL_KIND * PI / SLz

    nx0 = nx_global
    ny0 = ny_global
    nz0 = nz_global

    ! Set grid dimensions based on decomposition type
    select case (If_3dfft_decomp)
     case (0)  ! Block decomposition
       kx = nx_global / npx0
       ky = ny_global / npy0
       kz = nz_global / npz0
     case (1)  ! Slab decomposition
       kx = nx_global / npx0
       ky = ny_global / npz0
       kz = nz_global / npy0
     case (2)  ! Pencil decomposition
       kx = nx_global / npy0
       ky = ny_global / npz0
       kz = nz_global / npx0
    end select

    if (If_3dfft_decomp == 0) then
       ! Block decomposition: (x,y,z) -> (k_xi, k_eta, k_zeta)
       do k3 = 1, kz
          k_zeta = k3 + npz * nz - 1
          if (k_zeta >= nz0/2) k_zeta = k_zeta - nz0
          if (If_scheme == 0) kzeta_scaled = k_zeta * szz
          if (If_scheme == 1) kzeta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_zeta / nz0) / (SLz / nz0)

          do k2 = 1, ky
             k_eta = k2 + npy * ny - 1
             if (k_eta >= ny0/2) k_eta = k_eta - ny0
             if (If_scheme == 0) keta_scaled = k_eta * syy
             if (If_scheme == 1) keta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_eta / ny0) / (SLy / ny0)

             do k1 = 1, kx
                k_xi = k1 + npx * nx - 1
                if (k_xi >= nx0/2) cycle
                if (If_scheme == 0) kxi_scaled = k_xi * sxx
                if (If_scheme == 1) kxi_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_xi / nx0) / (SLx / nx0)

                k_squared_total = kxi_scaled**2 + keta_scaled**2 + kzeta_scaled**2

                if (k_squared_total > 0.0_ISFFT_REAL_KIND) then
                   inv_k_squared = -1.0_ISFFT_REAL_KIND / k_squared_total
                   phi_re_inout(k1, k2, k3) = phi_re_inout(k1, k2, k3) * inv_k_squared
                   phi_im_inout(k1, k2, k3) = phi_im_inout(k1, k2, k3) * inv_k_squared
                else
                   phi_re_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
                   phi_im_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
                endif
             enddo
          enddo
       enddo

    else
       ! Slab/Pencil decomposition: (z,x,y) -> (k_zeta, k_xi, k_eta)
       do k3 = 1, ky
          k_eta = k3 + npz * ky - 1
          if (k_eta >= ny0/2) k_eta = k_eta - ny0
          if (If_scheme == 0) keta_scaled = k_eta * syy
          if (If_scheme == 1) keta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_eta / ny0) / (SLy / ny0)

          do k2 = 1, kx
             k_xi = k2 + npy * kx - 1
             if (k_xi >= nx0/2) cycle
             if (If_scheme == 0) kxi_scaled = k_xi * sxx
             if (If_scheme == 1) kxi_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_xi / nx0) / (SLx / nx0)

             do k1 = 1, kz
                k_zeta = k1 + npx * kz - 1
                if (k_zeta >= nz0/2) k_zeta = k_zeta - nz0
                if (If_scheme == 0) kzeta_scaled = k_zeta * szz
                if (If_scheme == 1) kzeta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_zeta / nz0) / (SLz / nz0)

                k_squared_total = kxi_scaled**2 + keta_scaled**2 + kzeta_scaled**2

                if (k_squared_total > 0.0_ISFFT_REAL_KIND) then
                   inv_k_squared = -1.0_ISFFT_REAL_KIND / k_squared_total
                   phi_re_inout(k1, k2, k3) = phi_re_inout(k1, k2, k3) * inv_k_squared
                   phi_im_inout(k1, k2, k3) = phi_im_inout(k1, k2, k3) * inv_k_squared
                else
                   phi_re_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
                   phi_im_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
                endif
             enddo
          enddo
       enddo
    endif

 end subroutine Algebra_Periodic_Poisson_r2c
 !----------------------------------------------------------------------
 subroutine Algebra_Periodic_Poisson_r2c_serial(phi_re_inout, phi_im_inout)

    integer :: i, j, k, ierr
    integer :: kx, ky, kz, k1, k2, k3, k_xi, k_eta, k_zeta, nx0, ny0, nz0

    real(kind=ISFFT_REAL_KIND), intent(inout) :: phi_re_inout(:,:,:), phi_im_inout(:,:,:)
    real(kind=ISFFT_REAL_KIND) :: sxx, syy, szz
    real(kind=ISFFT_REAL_KIND) :: kxi_scaled, keta_scaled, kzeta_scaled, k_squared_total, inv_k_squared
    real(kind=ISFFT_REAL_KIND), parameter :: PI = 3.1415926535897932384626433832795_ISFFT_REAL_KIND

    ! Pre-compute constants
    sxx = 2.0_ISFFT_REAL_KIND * PI / SLx
    syy = 2.0_ISFFT_REAL_KIND * PI / SLy
    szz = 2.0_ISFFT_REAL_KIND * PI / SLz

    nx0 = nx_global
    ny0 = ny_global
    nz0 = nz_global

    kx = nx_global / npx0
    ky = ny_global / npy0
    kz = nz_global / npz0


    do k3 = 1, kz
       k_zeta = k3 + npz * nz - 1
       if (k_zeta >= nz0/2) k_zeta = k_zeta - nz0
       if (If_scheme == 0) kzeta_scaled = k_zeta * szz
       if (If_scheme == 1) kzeta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_zeta / nz0) / (SLz / nz0)

       do k2 = 1, ky
          k_eta = k2 + npy * ny - 1
          if (k_eta >= ny0/2) k_eta = k_eta - ny0
          if (If_scheme == 0) keta_scaled = k_eta * syy
          if (If_scheme == 1) keta_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_eta / ny0) / (SLy / ny0)

          do k1 = 1, kx
             k_xi = k1 + npx * nx - 1
             if (k_xi >= nx0/2) cycle
             if (If_scheme == 0) kxi_scaled = k_xi * sxx
             if (If_scheme == 1) kxi_scaled = 2.0_ISFFT_REAL_KIND * sin(PI * k_xi / nx0) / (SLx / nx0)

             k_squared_total = kxi_scaled**2 + keta_scaled**2 + kzeta_scaled**2

             if (k_squared_total > 0.0_ISFFT_REAL_KIND) then
                inv_k_squared = -1.0_ISFFT_REAL_KIND / k_squared_total
                phi_re_inout(k1, k2, k3) = phi_re_inout(k1, k2, k3) * inv_k_squared
                phi_im_inout(k1, k2, k3) = phi_im_inout(k1, k2, k3) * inv_k_squared
             else
                phi_re_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
                phi_im_inout(k1, k2, k3) = 0.0_ISFFT_REAL_KIND
             endif

          enddo
       enddo
    enddo

 end subroutine Algebra_Periodic_Poisson_r2c_serial
 !----------------------------------------------------------------------



end module ISFFT_Poisson_Algebra