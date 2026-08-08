clear; close all; clc;

% ============================================================
% Post-processing for 3D periodic Poisson tests
%
% If_case:
%   0: trigonometric periodic case
%   1: polynomial periodic case
% ============================================================

If_case = 0;   % 0: trigonometric periodic case
               % 1: polynomial periodic case

data = importdata('Poisson_3D_IO.dat');

Lx = 4; Ly = 4; Lz = 4;
nx = 128; ny = 128; nz = 128;

dx = Lx / nx;
dy = Ly / ny;
dz = Lz / nz;
dV = dx * dy * dz;

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
z = 0:dz:Lz-dz;

pp = reshape(data(:), [nx, ny, nz]);

switch If_case
    case 0
        % phi = sin(2*pi*z/Lz) * (sin(2*pi*x/Lx) + cos(4*pi*y/Ly))
        sx = sin(2*pi*x(:)/Lx);
        cy = cos(4*pi*y(:)/Ly);
        sz = sin(2*pi*z(:)/Lz);

        A = sx * sz.';
        B = cy * sz.';
        u = permute(A, [1, 3, 2]) + permute(B, [3, 1, 2]);

    case 1
        % phi = x^2(Lx-x)^2 y^2(Ly-y)^2 z^2(Lz-z)^2
        px = x(:).^2 .* (Lx - x(:)).^2;
        py = y(:).^2 .* (Ly - y(:)).^2;
        pz = z(:).^2 .* (Lz - z(:)).^2;

        Px = reshape(px, [nx, 1, 1]);
        Py = reshape(py, [1, ny, 1]);
        Pz = reshape(pz, [1, 1, nz]);

        u = Px .* Py .* Pz;
        u = u - mean(u(:));

    otherwise
        error('Unknown If_case. Please choose 0 or 1.');
end

pp0 = pp - mean(pp(:));
err = abs(pp0 - u);

Error_L2 = sqrt(sum(err(:).^2) * dV);
Error_Linf = max(err(:));

fprintf('Error_L2 = %.16e\n', Error_L2);
fprintf('Error_Linf = %.16e\n', Error_Linf);

ix = round(nx / 2);
iy = round(ny / 2);

figure;
plot(z, squeeze(u(ix,iy,:)), 'k-', 'LineWidth', 2); hold on;
plot(z, squeeze(pp0(ix,iy,:)), 'o', 'LineWidth', 1.2, 'MarkerSize', 5, ...
    'Color', [0.7608, 0.5059, 0.4863]);

xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\phi(x_0,y_0,z)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Theory', '$ISFFT$', 'Location', 'northeast', ...
    'Orientation', 'vertical', 'Box', 'off', ...
    'FontSize', 11, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
grid on;
set(gcf, 'Position', [100, 100, 560, 420]);
set(gca, 'LineWidth', 1.3, 'FontSize', 12);
