clear; close all; clc;

% ============================================================
% Post-processing for 3D Dirichlet Poisson tests
% Outputs:
%   Smooth cases:
%     1) Error_L2
%     2) Error_Linf
%   Dirac case:
%     no smooth L2/Linf error is printed
% and plots one center-line/source-line profile.
%
% Governing equation:
%   nabla^2 phi = f
%
% Boundary condition:
%   homogeneous Dirichlet for cases 0-5 and case 7
%   near-zero-at-boundary profile window for case 6
%
% Case 7:
%   3D discrete Dirac delta source under homogeneous Dirichlet BC.
%   Source point follows the Fortran setting:
%       isrc = nx_global / 2;
%       jsrc = ny_global / 2;
%       ksrc = nz_global / 2;
%       isrc = max(2, min(nx_global-1, isrc)); etc.
%   The source amplitude is:
%       f(isrc,jsrc,ksrc) = 1/(hx*hy*hz).
% ============================================================

If_case = 0;   % 0: trigonometric case 1
               % 1: trigonometric case 2
               % 2: quadratic polynomial case
               % 3: quartic polynomial case
               % 4: sextic polynomial case
               % 5: Gaussian-sine squared product case
               % 6: tanh-window product case (not exactly zero on boundary)
               % 7: 3D Dirac delta source, homogeneous Dirichlet BC

% Reference for case 7:
%   'discrete'   : exact solution of the second-order finite-difference
%                  Dirichlet problem with discrete sine eigenvalues.
%                  This is recommended if your Dirichlet solver uses the
%                  standard discrete Laplacian.
%   'continuous' : truncated continuous Green-function series with
%                  continuous eigenvalues. This is useful as a mathematical
%                  reference, but it does not exactly match the discrete
%                  operator near the singular source.
case7_reference_type = 'discrete';


% ---------------------------
% Basic parameters
% ---------------------------
data = importdata('Poisson_3D_IO.dat');

Lx = 4; Ly = 4; Lz = 4;
nx = 128; ny = 128; nz = 128;

case7_ref_mode_factor = 2;
case7_ref_max_modes = 384;

case7_continuous_modes_x = min(case7_ref_max_modes, case7_ref_mode_factor * (nx - 2));
case7_continuous_modes_y = min(case7_ref_max_modes, case7_ref_mode_factor * (ny - 2));
case7_continuous_modes_z = min(case7_ref_max_modes, case7_ref_mode_factor * (nz - 2));

case7_ref_line_n = max(1000, 4*nz);

% Dirichlet grid includes boundary points
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
z = linspace(0, Lz, nz);

hx = Lx / (nx - 1);
hy = Ly / (ny - 1);
hz = Lz / (nz - 1);

[X, Y, Z] = ndgrid(x, y, z);

% Numerical solution from solver output
pp = reshape(data(:), [nx, ny, nz]);

% ---------------------------
% Exact/reference solution for each case
% ---------------------------
src_info = struct();

switch If_case
    case 0
        % phi = sin(pi*x/Lx) * sin(2*pi*y/Ly) * sin(3*pi*z/Lz)
        u = sin(pi*X/Lx) .* sin(2*pi*Y/Ly) .* sin(3*pi*Z/Lz);
        case_name = 'Trig case 1';

    case 1
        % phi = sin(pi*x) * sin(pi*y) * sin(pi*z)
        u = sin(pi*X) .* sin(pi*Y) .* sin(pi*Z);
        case_name = 'Trig case 2';

    case 2
        % phi = x(Lx-x) y(Ly-y) z(Lz-z)
        u = X .* (Lx - X) .* Y .* (Ly - Y) .* Z .* (Lz - Z);
        case_name = 'Quadratic polynomial case';

    case 3
        % phi = [x^2 (Lx-x)^2] [y^2 (Ly-y)^2] [z^2 (Lz-z)^2]
        u = (X.^2) .* ((Lx - X).^2) .* ...
            (Y.^2) .* ((Ly - Y).^2) .* ...
            (Z.^2) .* ((Lz - Z).^2);
        case_name = 'Quartic polynomial case';

    case 4
        % phi = [x^3 (Lx-x)^3] [y^3 (Ly-y)^3] [z^3 (Lz-z)^3]
        u = (X.^3) .* ((Lx - X).^3) .* ...
            (Y.^3) .* ((Ly - Y).^3) .* ...
            (Z.^3) .* ((Lz - Z).^3);
        case_name = 'Sextic polynomial case';

    case 5
        % phi = prod_d exp(-((d-Ld/2)^2)/(2*sigma^2)) * sin(pi*d/Ld)^2
        sigma = 0.4;
        u = exp(-((X - 0.5*Lx).^2) / (2*sigma^2)) .* (sin(pi*X/Lx).^2) .* ...
            exp(-((Y - 0.5*Ly).^2) / (2*sigma^2)) .* (sin(pi*Y/Ly).^2) .* ...
            exp(-((Z - 0.5*Lz).^2) / (2*sigma^2)) .* (sin(pi*Z/Lz).^2);
        case_name = 'Gaussian-sine squared product case';

    case 6
        % phi = gx(x) * gy(y) * gz(z)
        % gx(x) = 0.5*(1+tanh((x-xstart)/s)) * 0.5*(1-tanh((x-xend)/s))
        % and similarly for y,z.
        % Note: this case is not exactly zero on the boundary in general.
        s = 0.2;
        xstart = 0.25 * Lx; xend = 0.75 * Lx;
        ystart = 0.25 * Ly; yend = 0.75 * Ly;
        zstart = 0.25 * Lz; zend = 0.75 * Lz;

        gx = 0.5 * (1 + tanh((X - xstart) / s)) .* 0.5 .* (1 - tanh((X - xend) / s));
        gy = 0.5 * (1 + tanh((Y - ystart) / s)) .* 0.5 .* (1 - tanh((Y - yend) / s));
        gz = 0.5 * (1 + tanh((Z - zstart) / s)) .* 0.5 .* (1 - tanh((Z - zend) / s));

        u = gx .* gy .* gz;
        case_name = 'Tanh-window product case';

    case 7
        % ------------------------------------------------------------
        % 3D Dirac delta source under homogeneous Dirichlet BC
        %
        % Continuous problem:
        %   nabla^2 phi = delta(x-xi)delta(y-eta)delta(z-zeta),
        %   phi = 0 on boundary.
        %
        % Discrete source:
        %   f(isrc,jsrc,ksrc) = 1/(hx*hy*hz), others 0.
        % ------------------------------------------------------------
        isrc = floor(nx / 2);
        jsrc = floor(ny / 2);
        ksrc = floor(nz / 2);

        isrc = max(2, min(nx-1, isrc));
        jsrc = max(2, min(ny-1, jsrc));
        ksrc = max(2, min(nz-1, ksrc));

        xi   = x(isrc);
        eta  = y(jsrc);
        zeta = z(ksrc);
        dV = hx * hy * hz;
        fsrc = 1.0 / dV;

        src_info.isrc = isrc;
        src_info.jsrc = jsrc;
        src_info.ksrc = ksrc;
        src_info.xi = xi;
        src_info.eta = eta;
        src_info.zeta = zeta;
        src_info.dV = dV;
        src_info.fsrc = fsrc;

        switch lower(case7_reference_type)
            case 'discrete'
                u = dirichlet_delta_green_discrete(nx, ny, nz, Lx, Ly, Lz, isrc, jsrc, ksrc);
                case_name = '3D discrete Dirac source, discrete Green reference';

            case 'continuous'
                u = dirichlet_delta_green_continuous(x, y, z, Lx, Ly, Lz, xi, eta, zeta, ...
                    case7_continuous_modes_x, case7_continuous_modes_y, case7_continuous_modes_z);
                case_name = sprintf('3D Dirac source, continuous Green reference, Mx=%d, My=%d, Mz=%d', ...
                    case7_continuous_modes_x, case7_continuous_modes_y, case7_continuous_modes_z);

            otherwise
                error('Unknown case7_reference_type. Use ''discrete'' or ''continuous''.');
        end

    otherwise
        error('Unknown If_case. Please choose 0, 1, 2, 3, 4, 5, 6, or 7.');
end

% ---------------------------
% Errors
% ---------------------------
err = abs(pp - u);

MAE = mean(err(:));
MaxError = max(err(:));
if If_case ~= 7
    dV = hx * hy * hz;
    Error_L2 = sqrt(sum(err(:).^2) * dV);
    Error_Linf = MaxError;
end

MaxBoundaryError = max([ ...
    max(abs(pp(1,:,:)   - u(1,:,:)),   [], 'all'), ...
    max(abs(pp(nx,:,:)  - u(nx,:,:)),  [], 'all'), ...
    max(abs(pp(:,1,:)   - u(:,1,:)),   [], 'all'), ...
    max(abs(pp(:,ny,:)  - u(:,ny,:)),  [], 'all'), ...
    max(abs(pp(:,:,1)   - u(:,:,1)),   [], 'all'), ...
    max(abs(pp(:,:,nz)  - u(:,:,nz)),  [], 'all')  ...
    ]);

% For the Dirac case, the point source is singular in the continuous problem.
% Even for the discrete reference, it is useful to report errors away from the
% source neighborhood.
if If_case == 7
    [II, JJ, KK] = ndgrid(1:nx, 1:ny, 1:nz);
    rr = sqrt(((II-src_info.isrc)*hx).^2 + ((JJ-src_info.jsrc)*hy).^2 + ((KK-src_info.ksrc)*hz).^2);
    mask_far = rr >= 3.0 * max([hx, hy, hz]);
    MAE_far = mean(err(mask_far));
    MaxError_far = max(err(mask_far));

    f_total = src_info.fsrc * src_info.dV;
    flux_num = boundary_flux_box(pp, x, y, z);
    flux_ref = boundary_flux_box(u, x, y, z);

    max_phi = max(pp(:));
    min_phi = min(pp(:));
    [~, min_id] = min(pp(:));
    [imin_phi, jmin_phi, kmin_phi] = ind2sub(size(pp), min_id);

    q = -pp;
    mono_tol = 1.0e-10 * max(abs(q(:)));
    mono_xm = monotonic_violations(q(src_info.isrc:-1:1, src_info.jsrc, src_info.ksrc), mono_tol);
    mono_xp = monotonic_violations(q(src_info.isrc:nx, src_info.jsrc, src_info.ksrc), mono_tol);
    mono_ym = monotonic_violations(squeeze(q(src_info.isrc, src_info.jsrc:-1:1, src_info.ksrc)), mono_tol);
    mono_yp = monotonic_violations(squeeze(q(src_info.isrc, src_info.jsrc:ny, src_info.ksrc)), mono_tol);
    mono_zm = monotonic_violations(squeeze(q(src_info.isrc, src_info.jsrc, src_info.ksrc:-1:1)), mono_tol);
    mono_zp = monotonic_violations(squeeze(q(src_info.isrc, src_info.jsrc, src_info.ksrc:nz)), mono_tol);

    r_boundary = min([src_info.xi, Lx-src_info.xi, src_info.eta, Ly-src_info.eta, src_info.zeta, Lz-src_info.zeta]);
    rmin_profile = 2.0 * max([hx, hy, hz]);
    rmax_profile = 0.9 * r_boundary;
    if rmax_profile > rmin_profile
        [rbin_num, prof_num] = radial_shell_profile(rr, -pp, rmin_profile, rmax_profile, 48);
        [rbin_ref, prof_ref] = radial_shell_profile(rr, -u,  rmin_profile, rmax_profile, 48);
        prof_free = 1.0 ./ (4.0*pi*rbin_num);

        fit_mask = isfinite(rbin_num) & isfinite(prof_num) & (prof_num > 0) & ...
                   (rbin_num >= 4.0*max([hx, hy, hz])) & (rbin_num <= 0.45*r_boundary);
        if nnz(fit_mask) >= 3
            pfit = polyfit(log(rbin_num(fit_mask)), log(prof_num(fit_mask)), 1);
            radial_loglog_slope = pfit(1);
        else
            radial_loglog_slope = NaN;
        end
    else
        rbin_num = []; prof_num = []; rbin_ref = []; prof_ref = []; prof_free = [];
        radial_loglog_slope = NaN;
    end
end

% ---------------------------
% Output
% ---------------------------
if If_case ~= 7
    fprintf('Error_L2 = %.16e\n', Error_L2);
    fprintf('Error_Linf = %.16e\n', Error_Linf);
end

% ---------------------------
% Profile plot
% ---------------------------
if If_case == 7
    ix = src_info.isrc;
    iy = src_info.jsrc;
else
    ix = round((nx + 1) / 2);
    iy = round((ny + 1) / 2);
end

figure;
plot(z, squeeze(u(ix,iy,:)), 'k-', 'LineWidth', 2); hold on;
plot(z, squeeze(pp(ix,iy,:)), 'o', 'LineWidth', 1.2, 'MarkerSize', 5, ...
    'Color', [0.7608, 0.5059, 0.4863]);

if If_case == 7
    xline(src_info.zeta, '--', '$z=\zeta$', 'Interpreter', 'latex', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
end

xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\phi(x_0,y_0,z)$', 'Interpreter', 'latex', 'FontSize', 16);
% title(sprintf('%s, If\\_case = %d', case_name, If_case), 'FontSize', 13);
legend('Reference', '$ISFFT$', 'Location', 'best', ...
    'Orientation', 'vertical', 'Box', 'off', ...
    'FontSize', 11, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

grid on;
set(gcf, 'Position', [100, 100, 560, 420]);
set(gca, 'LineWidth', 1.3, 'FontSize', 12);

% ---------------------------
% Optional slice plot for case 7
% ---------------------------
if If_case == 7
    figure;
    iz = src_info.ksrc;
    imagesc(x, y, squeeze(pp(:,:,iz))'); axis image; set(gca, 'YDir', 'normal');
    hold on;
    plot(src_info.xi, src_info.eta, 'wo', 'MarkerSize', 7, 'LineWidth', 1.5);
    colorbar;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
    title(sprintf('Numerical solution slice at z = %.4f', src_info.zeta), 'FontSize', 13);
    set(gcf, 'Position', [700, 100, 520, 430]);
    set(gca, 'LineWidth', 1.3, 'FontSize', 12);

    if ~isempty(rbin_num)
        figure;
        loglog(rbin_num, prof_num, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5); hold on;
        loglog(rbin_ref, prof_ref, 'k-', 'LineWidth', 1.8);
        loglog(rbin_num, prof_free, '--', 'LineWidth', 1.5);
        xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 16);
        ylabel('$-\langle \phi \rangle_{r}$', 'Interpreter', 'latex', 'FontSize', 16);
        title(sprintf('Radial decay check, slope = %.3f', radial_loglog_slope), 'FontSize', 13);
        legend('$ISFFT$', 'Reference', '$1/(4\pi r)$', 'Location', 'best', ...
            'Box', 'off', 'Interpreter', 'latex', 'FontSize', 11);
        grid on;
        set(gcf, 'Position', [700, 560, 540, 430]);
        set(gca, 'LineWidth', 1.3, 'FontSize', 12);
    end
end

% Optional:
% print(sprintf('profile_case_%d', If_case), '-dpng', '-r300');

% ============================================================
% Local functions
% ============================================================
function flux = boundary_flux_box(phi, x, y, z)
    hx = x(2) - x(1);
    hy = y(2) - y(1);
    hz = z(2) - z(1);
    nx = numel(x); ny = numel(y); nz = numel(z);

    dphidn_x0 = squeeze(-(phi(2,:,:) - phi(1,:,:)) / hx);
    dphidn_x1 = squeeze( (phi(nx,:,:) - phi(nx-1,:,:)) / hx);
    dphidn_y0 = squeeze(-(phi(:,2,:) - phi(:,1,:)) / hy);
    dphidn_y1 = squeeze( (phi(:,ny,:) - phi(:,ny-1,:)) / hy);
    dphidn_z0 = squeeze(-(phi(:,:,2) - phi(:,:,1)) / hz);
    dphidn_z1 = squeeze( (phi(:,:,nz) - phi(:,:,nz-1)) / hz);

    flux_x0 = trapz(y, trapz(z, dphidn_x0, 2));
    flux_x1 = trapz(y, trapz(z, dphidn_x1, 2));
    flux_y0 = trapz(x, trapz(z, dphidn_y0, 2));
    flux_y1 = trapz(x, trapz(z, dphidn_y1, 2));
    flux_z0 = trapz(x, trapz(y, dphidn_z0, 2));
    flux_z1 = trapz(x, trapz(y, dphidn_z1, 2));

    flux = flux_x0 + flux_x1 + flux_y0 + flux_y1 + flux_z0 + flux_z1;
end

function nbad = monotonic_violations(v, tol)
    v = v(:);
    if numel(v) <= 1
        nbad = 0;
        return;
    end
    nbad = nnz(diff(v) > tol);
end

function [rbin, prof] = radial_shell_profile(rr, val, rmin, rmax, nbins)
    edges = logspace(log10(rmin), log10(rmax), nbins+1);
    rbin = sqrt(edges(1:end-1) .* edges(2:end));
    prof = NaN(size(rbin));
    for ib = 1:nbins
        mask = (rr >= edges(ib)) & (rr < edges(ib+1));
        vv = val(mask);
        vv = vv(isfinite(vv) & (vv > 0));
        if ~isempty(vv)
            prof(ib) = mean(vv);
        end
    end
end

function u = dirichlet_delta_green_discrete(nx, ny, nz, Lx, Ly, Lz, isrc, jsrc, ksrc)
    % Exact discrete Green function for the standard second-order finite-
    % difference Laplacian with homogeneous Dirichlet boundary conditions.
    %
    % Grid includes boundary points. Interior indices are 2:nx-1, etc.
    % Source follows:
    %   f(isrc,jsrc,ksrc) = 1/(hx*hy*hz).

    hx = Lx / (nx - 1);
    hy = Ly / (ny - 1);
    hz = Lz / (nz - 1);
    dV = hx * hy * hz;
    fsrc = 1.0 / dV;

    nix = nx - 2;
    niy = ny - 2;
    niz = nz - 2;

    ax_src = isrc - 1;   % interior index, 1:nix
    ay_src = jsrc - 1;
    az_src = ksrc - 1;

    p = (1:nix)';
    q = (1:niy)';
    r = (1:niz)';

    % Interior physical/grid indices a=1:nix correspond to global i=a+1.
    ax = (1:nix)';
    ay = (1:niy)';
    az = (1:niz)';

    Sx = sin(pi * ax * p' / (nx - 1));
    Sy = sin(pi * ay * q' / (ny - 1));
    Sz = sin(pi * az * r' / (nz - 1));

    sx_src = sin(pi * p * ax_src / (nx - 1));
    sy_src = sin(pi * q * ay_src / (ny - 1));
    sz_src = sin(pi * r * az_src / (nz - 1));

    lamx = (2.0*cos(pi*p/(nx-1)) - 2.0) / hx^2;
    lamy = (2.0*cos(pi*q/(ny-1)) - 2.0) / hy^2;
    lamz = (2.0*cos(pi*r/(nz-1)) - 2.0) / hz^2;

    coeff = zeros(nix, niy, niz);
    norm_factor = 8.0 / ((nx - 1) * (ny - 1) * (nz - 1));

    for kk = 1:niz
        for jj = 1:niy
            denom = lamx + lamy(jj) + lamz(kk);
            coeff(:,jj,kk) = norm_factor * fsrc * sx_src .* sy_src(jj) .* sz_src(kk) ./ denom;
        end
    end

    uint = sine_tensor_reconstruct(coeff, Sx, Sy, Sz);

    u = zeros(nx, ny, nz);
    u(2:nx-1, 2:ny-1, 2:nz-1) = uint;
end

function u = dirichlet_delta_green_continuous(x, y, z, Lx, Ly, Lz, xi, eta, zeta, Mx, My, Mz)
    % Truncated continuous Dirichlet Green-function series for
    %   nabla^2 u = delta(x-xi)delta(y-eta)delta(z-zeta), u=0 on boundary.
    % The sign is negative because the governing equation is nabla^2 u = f.

    nx = numel(x); ny = numel(y); nz = numel(z);
    m = (1:Mx)';
    n = (1:My)';
    l = (1:Mz)';

    Sx = sin(pi * x(:) * m' / Lx);
    Sy = sin(pi * y(:) * n' / Ly);
    Sz = sin(pi * z(:) * l' / Lz);

    sx_src = sin(pi * m * xi / Lx);
    sy_src = sin(pi * n * eta / Ly);
    sz_src = sin(pi * l * zeta / Lz);

    lamx = (pi*m/Lx).^2;
    lamy = (pi*n/Ly).^2;
    lamz = (pi*l/Lz).^2;

    coeff = zeros(Mx, My, Mz);
    factor = -8.0 / (Lx * Ly * Lz);

    for kk = 1:Mz
        for jj = 1:My
            denom = lamx + lamy(jj) + lamz(kk);
            coeff(:,jj,kk) = factor * sx_src .* sy_src(jj) .* sz_src(kk) ./ denom;
        end
    end

    u = sine_tensor_reconstruct(coeff, Sx, Sy, Sz);
    u = reshape(u, [nx, ny, nz]);
end

function u = sine_tensor_reconstruct(coeff, Sx, Sy, Sz)
    % Reconstruct u(i,j,k) = sum_{p,q,r} coeff(p,q,r) Sx(i,p) Sy(j,q) Sz(k,r)
    % using three matrix multiplications.

    [mx, my, mz] = size(coeff);
    nx = size(Sx, 1);
    ny = size(Sy, 1);
    nz = size(Sz, 1);

    tmp = Sx * reshape(coeff, mx, my*mz);       % nx x (my*mz)
    tmp = reshape(tmp, nx, my, mz);             % nx x my x mz

    tmp = permute(tmp, [2, 1, 3]);              % my x nx x mz
    tmp = Sy * reshape(tmp, my, nx*mz);         % ny x (nx*mz)
    tmp = reshape(tmp, ny, nx, mz);             % ny x nx x mz
    tmp = permute(tmp, [2, 1, 3]);              % nx x ny x mz

    tmp = permute(tmp, [3, 1, 2]);              % mz x nx x ny
    tmp = Sz * reshape(tmp, mz, nx*ny);         % nz x (nx*ny)
    tmp = reshape(tmp, nz, nx, ny);             % nz x nx x ny
    u = permute(tmp, [2, 3, 1]);                % nx x ny x nz
end
