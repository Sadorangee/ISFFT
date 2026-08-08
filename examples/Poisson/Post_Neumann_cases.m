clear; close all; clc;

% ============================================================
% Post-processing for 3D homogeneous Neumann Poisson tests
%
% Governing equation:
%   nabla^2 phi = f
%
% Boundary condition:
%   dphi/dn = 0 on all six boundary faces
%
% Grid convention for homogeneous Neumann with DCT-I:
%   The physical grid includes both boundary points,
%       x_i = (i-1)*Lx/(nx-1), i = 1,...,nx.
%   This is different from a Dirichlet DST solve, where only interior
%   unknowns participate in the spectral algebra.
%
% If_case:
%   0: trigonometric smooth Neumann case
%   1: polynomial smooth Neumann case
%   2: discrete Dirac source with zero-mean correction
%
% For cases 0 and 1, the exact solution is known. Since a pure Neumann
% solution is determined only up to an additive constant, the numerical
% solution is aligned to the exact solution before computing L2/Linf.
%
% For case 2, no smooth classical exact solution is used. The diagnostics
% focus on source balance, residual away from the singular source, source
% location, and near-source radial decay.
% ============================================================

If_case = 0;   % 0: trigonometric Neumann case
               % 1: polynomial Neumann case
               % 2: discrete Dirac Neumann case

% ---------------------------
% Basic parameters
% ---------------------------
data = importdata('Poisson_3D_IO.dat');

Lx = 4; Ly = 4; Lz = 4;
nx = 256; ny = 256; nz = 256;

% Neumann DCT-I grid includes boundary points.
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
z = linspace(0, Lz, nz);

hx = Lx / (nx - 1);
hy = Ly / (ny - 1);
hz = Lz / (nz - 1);
hmax = max([hx, hy, hz]);

[X, Y, Z] = ndgrid(x, y, z);

% Numerical solution from solver output
pp = reshape(data(:), [nx, ny, nz]);

% ---------------------------
% Exact/reference solution and source
% ---------------------------
has_exact_solution = true;
src_info = struct();

switch If_case
    case 0
        % ------------------------------------------------------------
        % Trigonometric Neumann case
        %
        % phi = cos(pi*x) cos(2*pi*y) cos(3*pi*z)
        % f   = -[pi^2 + (2*pi)^2 + (3*pi)^2] phi
        %
        % This satisfies homogeneous Neumann BC when Lx, Ly, and Lz
        % are integers, e.g. Lx = Ly = Lz = 4 in this file.
        % ------------------------------------------------------------
        mx = 1; my = 2; mz = 3;

        u = cos(mx*pi*X) .* cos(my*pi*Y) .* cos(mz*pi*Z);

        lambda_cont = (mx*pi)^2 + (my*pi)^2 + (mz*pi)^2;
        f = -lambda_cont * u;

        case_name = 'Trigonometric Neumann case';

    case 1
        % ------------------------------------------------------------
        % Polynomial Neumann case
        %
        % phi = gx(x) gy(y) gz(z)
        % gx = x^2 (Lx-x)^2, etc.
        % This satisfies dphi/dn = 0 on all boundaries.
        % ------------------------------------------------------------
        gx = X.^2 .* (Lx - X).^2;
        gy = Y.^2 .* (Ly - Y).^2;
        gz = Z.^2 .* (Lz - Z).^2;

        g2x = 12*X.^2 - 12*Lx*X + 2*Lx^2;
        g2y = 12*Y.^2 - 12*Ly*Y + 2*Ly^2;
        g2z = 12*Z.^2 - 12*Lz*Z + 2*Lz^2;

        u = gx .* gy .* gz;
        f = g2x .* gy .* gz + gx .* g2y .* gz + gx .* gy .* g2z;

        case_name = 'Polynomial Neumann case';

    case 2
        % ------------------------------------------------------------
        % Discrete Dirac source for pure Neumann Poisson
        %
        % Continuous form:
        %   nabla^2 phi = delta(x-x0,y-y0,z-z0) - 1/|Omega|
        %
        % Discrete source used here:
        %   f = delta_h - mean(delta_h)
        %   delta_h = 1/dV at the source grid point.
        %
        % This case has a singular Green-function response, so no smooth
        % exact solution is used for MAE.
        % ------------------------------------------------------------
        has_exact_solution = false;

        dV = hx * hy * hz;

        isrc = floor(nx / 2);
        jsrc = floor(ny / 2);
        ksrc = floor(nz / 2);

        isrc = max(1, min(nx, isrc));
        jsrc = max(1, min(ny, jsrc));
        ksrc = max(1, min(nz, ksrc));

        f = -ones(nx, ny, nz) / (dV * nx * ny * nz);
        f(isrc, jsrc, ksrc) = f(isrc, jsrc, ksrc) + 1.0 / dV;

        u = [];

        src_info.isrc = isrc;
        src_info.jsrc = jsrc;
        src_info.ksrc = ksrc;
        src_info.xi = x(isrc);
        src_info.eta = y(jsrc);
        src_info.zeta = z(ksrc);
        src_info.dV = dV;

        case_name = 'Discrete Dirac Neumann case';

    otherwise
        error('Unknown If_case. Please choose 0, 1, or 2.');
end

% ---------------------------
% Pure Neumann constant alignment
% ---------------------------
if has_exact_solution
    constant_offset = mean(pp(:) - u(:));
    pp_aligned = pp - constant_offset;
else
    pp_aligned = pp - mean(pp(:));
end

% ---------------------------
% Residual diagnostic
% ---------------------------
% This residual is evaluated with a standard second-order Neumann
% finite-difference Laplacian. It is used only as a post-processing
% diagnostic.
lap_pp = neumann_laplacian_2nd(pp_aligned, hx, hy, hz);
residual = lap_pp - f;

% ---------------------------
% Output and plots
% ---------------------------
if has_exact_solution
    err = abs(pp_aligned - u);

    dV = hx * hy * hz;
    Error_L2 = sqrt(sum(err(:).^2) * dV);
    Error_Linf = max(err(:));
    MAE = mean(err(:));
    Residual_MAE = mean(abs(residual(:)));

    fprintf('Error_L2 = %.16e\n', Error_L2);
    fprintf('Error_Linf = %.16e\n', Error_Linf);

    ix = round((nx + 1) / 2);
    iy = round((ny + 1) / 2);

    figure;
    plot(z, squeeze(u(ix,iy,:)), 'k-', 'LineWidth', 2); hold on;
    plot(z, squeeze(pp_aligned(ix,iy,:)), 'o', 'LineWidth', 1.2, 'MarkerSize', 5, ...
        'Color', [0.7608, 0.5059, 0.4863]);

    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\phi(x_0,y_0,z)$', 'Interpreter', 'latex', 'FontSize', 16);
    legend('Exact', '$ISFFT$', 'Location', 'best', ...
        'Orientation', 'vertical', 'Box', 'off', ...
        'FontSize', 11, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    grid on;
    set(gcf, 'Position', [100, 100, 560, 420]);
    set(gca, 'LineWidth', 1.3, 'FontSize', 12);

    iz = round((nz + 1) / 2);
    figure;
    imagesc(x, y, squeeze(err(:,:,iz))'); axis image; set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
%     title(sprintf('Absolute error slice at z = %.4f', z(iz)), 'FontSize', 13);
    set(gcf, 'Position', [700, 100, 520, 430]);
    set(gca, 'LineWidth', 1.3, 'FontSize', 12);

else
    [II, JJ, KK] = ndgrid(1:nx, 1:ny, 1:nz);
    rr = sqrt(((II-src_info.isrc)*hx).^2 + ...
              ((JJ-src_info.jsrc)*hy).^2 + ...
              ((KK-src_info.ksrc)*hz).^2);

    mask_far = rr >= 3.0 * hmax;
    Residual_MAE_far = mean(abs(residual(mask_far)));

    SourceBalance = sum(f(:)) * src_info.dV;

    [min_phi, min_id] = min(pp_aligned(:));
    [imin_phi, jmin_phi, kmin_phi] = ind2sub(size(pp_aligned), min_id);
    SourceIsMinimum = (imin_phi == src_info.isrc) && ...
                      (jmin_phi == src_info.jsrc) && ...
                      (kmin_phi == src_info.ksrc);

    r_boundary = min([src_info.xi, Lx-src_info.xi, ...
                      src_info.eta, Ly-src_info.eta, ...
                      src_info.zeta, Lz-src_info.zeta]);
    rmin_profile = 2.0 * hmax;
    rmax_profile = 0.9 * r_boundary;

    if rmax_profile > rmin_profile
        [rbin_num, prof_num] = radial_shell_profile(rr, -pp_aligned, rmin_profile, rmax_profile, 48);
        prof_free = 1.0 ./ (4.0*pi*rbin_num);

        fit_mask = isfinite(rbin_num) & isfinite(prof_num) & (prof_num > 0) & ...
                   (rbin_num >= 4.0*hmax) & (rbin_num <= 0.45*r_boundary);
        if nnz(fit_mask) >= 3
            pfit = polyfit(log(rbin_num(fit_mask)), log(prof_num(fit_mask)), 1);
            RadialLogLogSlope = pfit(1);
        else
            RadialLogLogSlope = NaN;
        end
    else
        rbin_num = [];
        prof_num = [];
        prof_free = [];
        RadialLogLogSlope = NaN;
    end

    ix = src_info.isrc;
    iy = src_info.jsrc;

    figure;
    plot(z, squeeze(pp_aligned(ix,iy,:)), 'o-', 'LineWidth', 1.2, 'MarkerSize', 5, ...
        'Color', [0.7608, 0.5059, 0.4863]); hold on;
    xline(src_info.zeta, '--', '$z=\zeta$', 'Interpreter', 'latex', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\phi(x_s,y_s,z)$', 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    set(gcf, 'Position', [100, 100, 560, 420]);
    set(gca, 'LineWidth', 1.3, 'FontSize', 12);

    figure;
    iz = src_info.ksrc;
    imagesc(x, y, squeeze(pp_aligned(:,:,iz))'); axis image; set(gca, 'YDir', 'normal');
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
        loglog(rbin_num, prof_free, '--', 'LineWidth', 1.5);
        xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 16);
        ylabel('$-\langle \phi \rangle_r$', 'Interpreter', 'latex', 'FontSize', 16);
        title(sprintf('Radial decay check, slope = %.3f', RadialLogLogSlope), 'FontSize', 13);
        legend('$ISFFT$', '$1/(4\pi r)$', 'Location', 'best', ...
            'Box', 'off', 'Interpreter', 'latex', 'FontSize', 11);
        grid on;
        set(gcf, 'Position', [700, 560, 540, 430]);
        set(gca, 'LineWidth', 1.3, 'FontSize', 12);
    end
end

% Optional:
% print(sprintf('neumann_case_%d', If_case), '-dpng', '-r300');

% ============================================================
% Local functions
% ============================================================
function lap = neumann_laplacian_2nd(phi, hx, hy, hz)
    [nx, ny, nz] = size(phi);
    lap = zeros(nx, ny, nz);

    % x-direction
    lap(2:nx-1,:,:) = lap(2:nx-1,:,:) + ...
        (phi(3:nx,:,:) - 2*phi(2:nx-1,:,:) + phi(1:nx-2,:,:)) / hx^2;
    lap(1,:,:) = lap(1,:,:) + ...
        2*(phi(2,:,:) - phi(1,:,:)) / hx^2;
    lap(nx,:,:) = lap(nx,:,:) + ...
        2*(phi(nx-1,:,:) - phi(nx,:,:)) / hx^2;

    % y-direction
    lap(:,2:ny-1,:) = lap(:,2:ny-1,:) + ...
        (phi(:,3:ny,:) - 2*phi(:,2:ny-1,:) + phi(:,1:ny-2,:)) / hy^2;
    lap(:,1,:) = lap(:,1,:) + ...
        2*(phi(:,2,:) - phi(:,1,:)) / hy^2;
    lap(:,ny,:) = lap(:,ny,:) + ...
        2*(phi(:,ny-1,:) - phi(:,ny,:)) / hy^2;

    % z-direction
    lap(:,:,2:nz-1) = lap(:,:,2:nz-1) + ...
        (phi(:,:,3:nz) - 2*phi(:,:,2:nz-1) + phi(:,:,1:nz-2)) / hz^2;
    lap(:,:,1) = lap(:,:,1) + ...
        2*(phi(:,:,2) - phi(:,:,1)) / hz^2;
    lap(:,:,nz) = lap(:,:,nz) + ...
        2*(phi(:,:,nz-1) - phi(:,:,nz)) / hz^2;
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
