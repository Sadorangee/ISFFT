clear; close all; clc;

% ========================================================================
% Unified post-processing for the smooth manufactured Poisson test cases.
%
% Match the settings below to Poisson_solver_UI.in before running this file.
% Boundary values: 0 = Periodic, 1 = Dirichlet, 2 = Neumann.
%
% If_case is used only for the 3D pure-boundary families:
%   Boundary = [0 0 0]: If_case = 0:1 (periodic tests)
%   Boundary = [1 1 1]: If_case = 0:6 (smooth Dirichlet tests)
%   Boundary = [2 2 2]: If_case = 0:1 (smooth Neumann tests)
% It is ignored for 1D, 2D, and 3D mixed-boundary tests.
%
% Singular cases use their dedicated diagnostics instead:
%   3D Dirichlet If_case = 7 -> Post_Dirichlet_cases.m
%   3D Neumann   If_case = 2 -> Post_Neumann_cases.m
% ========================================================================

Dimension = 3;
Boundary = [0, 0, 0];
If_case = 0;

% These values must match Poisson_solver_UI.in.
Grid_Size = [128, 128, 128];
Domain_Size = [4.0, 4.0, 4.0];

if ~ismember(Dimension, 1:3)
    error('Dimension must be 1, 2, or 3.');
end
if numel(Boundary) ~= 3 || any(~ismember(Boundary, 0:2))
    error('Boundary must contain three values chosen from 0, 1, and 2.');
end
if any(Grid_Size < 1) || any(Domain_Size <= 0)
    error('Grid_Size must be positive and Domain_Size must be greater than zero.');
end
if Dimension == 1 && any(Grid_Size(2:3) ~= 1)
    error('Dimension = 1 requires Grid_Size(2:3) = [1, 1].');
elseif Dimension == 2 && Grid_Size(3) ~= 1
    error('Dimension = 2 requires Grid_Size(3) = 1.');
end

if Dimension == 3 && all(Boundary == 1) && If_case == 7
    error('Use Post_Dirichlet_cases.m for the singular Dirichlet case 7.');
elseif Dimension == 3 && all(Boundary == 2) && If_case == 2
    error('Use Post_Neumann_cases.m for the singular Neumann case 2.');
end

nx = Grid_Size(1); ny = Grid_Size(2); nz = Grid_Size(3);
Lx = Domain_Size(1); Ly = Domain_Size(2); Lz = Domain_Size(3);

x = grid_coordinates(nx, Lx, Boundary(1));
y = grid_coordinates(ny, Ly, Boundary(2));
z = grid_coordinates(nz, Lz, Boundary(3));

switch Dimension
    case 1
        output_file = 'Poisson_1D_IO.dat';
        data = importdata(output_file);
        assert_output_size(data, nx, output_file);
        pp = reshape(data(:), [nx, 1]);
        u = exact_1d(x, Lx, Boundary(1));
        u = reshape(u, [nx, 1]);
        case_name = sprintf('1D %s case', boundary_name(Boundary(1)));

    case 2
        output_file = 'Poisson_2D_IO.dat';
        data = importdata(output_file);
        assert_output_size(data, nx*ny, output_file);
        pp = reshape(data(:), [nx, ny]);
        [X, Y] = ndgrid(x, y);
        u = exact_2d(X, Y, Lx, Ly, Boundary(1:2));
        case_name = sprintf('2D %s-%s case', ...
            boundary_name(Boundary(1)), boundary_name(Boundary(2)));

    case 3
        output_file = 'Poisson_3D_IO.dat';
        data = importdata(output_file);
        assert_output_size(data, nx*ny*nz, output_file);
        pp = reshape(data(:), [nx, ny, nz]);
        [X, Y, Z] = ndgrid(x, y, z);
        [u, case_name] = exact_3d(X, Y, Z, Lx, Ly, Lz, Boundary, If_case);
end

% Pure periodic, pure Neumann, and periodic/Neumann mixed problems have a
% constant null mode. Align the numerical solution before computing errors.
if ~any(Boundary(1:Dimension) == 1)
    pp = pp - mean(pp(:)) + mean(u(:));
end

err = pp - u;
reference_norm = norm(u(:), 2);
if reference_norm <= eps(max(1.0, norm(u(:), inf)))
    error('The reference solution has zero L2 norm; a relative error is undefined.');
end

Relative_L2_Error = norm(err(:), 2) / reference_norm;
Linf_Error = norm(err(:), inf);

fprintf('%s\n', case_name);
fprintf('Relative_L2_Error = %.16e\n', Relative_L2_Error);
fprintf('Linf_Error        = %.16e\n', Linf_Error);

% Compare one representative profile along the last active coordinate.
switch Dimension
    case 1
        profile_coordinate = x;
        profile_exact = u;
        profile_isfft = pp;
        coordinate_label = '$x$';
        solution_label = '$\phi(x)$';

    case 2
        ix = round((nx + 1) / 2);
        profile_coordinate = y;
        profile_exact = squeeze(u(ix,:));
        profile_isfft = squeeze(pp(ix,:));
        coordinate_label = '$y$';
        solution_label = '$\phi(x_0,y)$';

    case 3
        ix = round((nx + 1) / 2);
        iy = round((ny + 1) / 2);
        profile_coordinate = z;
        profile_exact = squeeze(u(ix,iy,:));
        profile_isfft = squeeze(pp(ix,iy,:));
        coordinate_label = '$z$';
        solution_label = '$\phi(x_0,y_0,z)$';
end

figure;
plot(profile_coordinate, profile_exact, 'k-', 'LineWidth', 2); hold on;
plot(profile_coordinate, profile_isfft, 'o', 'LineWidth', 1.2, ...
    'MarkerSize', 5, 'Color', [0.7608, 0.5059, 0.4863]);
xlabel(coordinate_label, 'Interpreter', 'latex', 'FontSize', 16);
ylabel(solution_label, 'Interpreter', 'latex', 'FontSize', 16);
legend('Exact', '$ISFFT$', 'Location', 'best', 'Orientation', 'vertical', ...
    'Box', 'off', 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'Interpreter', 'latex');
grid on;
set(gcf, 'Position', [100, 100, 560, 420]);
set(gca, 'LineWidth', 1.3, 'FontSize', 12);

% ========================================================================
% Local functions
% ========================================================================
function coordinate = grid_coordinates(n, length_value, boundary_type)
    if n == 1
        coordinate = 0;
    elseif boundary_type == 0
        coordinate = (0:n-1) * length_value / n;
    else
        coordinate = linspace(0, length_value, n);
    end
end

function assert_output_size(data, expected_size, filename)
    if numel(data) ~= expected_size
        error('%s contains %d values, but the selected settings require %d.', ...
            filename, numel(data), expected_size);
    end
end

function name = boundary_name(boundary_type)
    names = {'Periodic', 'Dirichlet', 'Neumann'};
    name = names{boundary_type + 1};
end

function u = exact_1d(x, Lx, bx)
    switch bx
        case 0
            u = sin(2*pi*x/Lx);
        case 1
            u = sin(pi*x/Lx);
        case 2
            u = cos(pi*x/Lx);
    end
end

function u = exact_2d(X, Y, Lx, Ly, boundary)
    bx = boundary(1); by = boundary(2);
    key = 10*bx + by;

    switch key
        case 0
            u = sin(2*pi*X/Lx) + cos(4*pi*Y/Ly);
        case 11
            u = sin(pi*X/Lx) .* sin(2*pi*Y/Ly);
        case 22
            u = cos(pi*X/Lx) .* cos(2*pi*Y/Ly);
        case 1
            u = cos(2*pi*X/Lx) .* sin(pi*Y/Ly);
        case 2
            u = cos(2*pi*X/Lx) .* cos(pi*Y/Ly);
        case 12
            u = sin(pi*X/Lx) .* cos(pi*Y/Ly);
        otherwise
            u = axis_factor(X, Lx, bx) .* axis_factor(Y, Ly, by);
    end
end

function [u, case_name] = exact_3d(X, Y, Z, Lx, Ly, Lz, boundary, if_case)
    if all(boundary == 0)
        switch if_case
            case 0
                u = sin(2*pi*Z/Lz) .* ...
                    (sin(2*pi*X/Lx) + cos(4*pi*Y/Ly));
                case_name = '3D periodic trigonometric case';
            case 1
                u = X.^2 .* (Lx-X).^2 .* Y.^2 .* (Ly-Y).^2 .* ...
                    Z.^2 .* (Lz-Z).^2;
                u = u - mean(u(:));
                case_name = '3D periodic polynomial case';
            otherwise
                error('For 3D pure periodic tests, If_case must be 0 or 1.');
        end
        return;
    end

    if all(boundary == 1)
        switch if_case
            case 0
                u = sin(pi*X/Lx) .* sin(2*pi*Y/Ly) .* sin(3*pi*Z/Lz);
                case_name = '3D Dirichlet trigonometric case 0';
            case 1
                u = sin(pi*X) .* sin(pi*Y) .* sin(pi*Z);
                case_name = '3D Dirichlet trigonometric case 1';
            case 2
                u = X .* (Lx-X) .* Y .* (Ly-Y) .* Z .* (Lz-Z);
                case_name = '3D Dirichlet quadratic polynomial case';
            case 3
                u = X.^2 .* (Lx-X).^2 .* Y.^2 .* (Ly-Y).^2 .* ...
                    Z.^2 .* (Lz-Z).^2;
                case_name = '3D Dirichlet quartic polynomial case';
            case 4
                u = X.^3 .* (Lx-X).^3 .* Y.^3 .* (Ly-Y).^3 .* ...
                    Z.^3 .* (Lz-Z).^3;
                case_name = '3D Dirichlet sextic polynomial case';
            case 5
                sigma = 0.4;
                u = exp(-(X-Lx/2).^2/(2*sigma^2)) .* sin(pi*X/Lx).^2 .* ...
                    exp(-(Y-Ly/2).^2/(2*sigma^2)) .* sin(pi*Y/Ly).^2 .* ...
                    exp(-(Z-Lz/2).^2/(2*sigma^2)) .* sin(pi*Z/Lz).^2;
                case_name = '3D Dirichlet Gaussian-sine-squared case';
            case 6
                width = 0.2;
                gx = tanh_window(X, 0.25*Lx, 0.75*Lx, width);
                gy = tanh_window(Y, 0.25*Ly, 0.75*Ly, width);
                gz = tanh_window(Z, 0.25*Lz, 0.75*Lz, width);
                u = gx .* gy .* gz;
                case_name = '3D Dirichlet tanh-window case';
            otherwise
                error('For smooth 3D Dirichlet tests, If_case must be between 0 and 6.');
        end
        return;
    end

    if all(boundary == 2)
        switch if_case
            case 0
                u = cos(pi*X) .* cos(2*pi*Y) .* cos(3*pi*Z);
                case_name = '3D Neumann trigonometric case';
            case 1
                u = X.^2 .* (Lx-X).^2 .* Y.^2 .* (Ly-Y).^2 .* ...
                    Z.^2 .* (Lz-Z).^2;
                case_name = '3D Neumann polynomial case';
            otherwise
                error('For smooth 3D Neumann tests, If_case must be 0 or 1.');
        end
        return;
    end

    u = axis_factor(X, Lx, boundary(1)) .* ...
        axis_factor(Y, Ly, boundary(2)) .* ...
        axis_factor(Z, Lz, boundary(3));
    case_name = sprintf('3D %s-%s-%s mixed-boundary case', ...
        boundary_name(boundary(1)), boundary_name(boundary(2)), ...
        boundary_name(boundary(3)));
end

function value = axis_factor(coordinate, length_value, boundary_type)
    switch boundary_type
        case 0
            value = cos(2*pi*coordinate/length_value);
        case 1
            value = sin(pi*coordinate/length_value);
        case 2
            value = cos(pi*coordinate/length_value);
    end
end

function value = tanh_window(coordinate, start_value, end_value, width)
    value = 0.5 .* (1 + tanh((coordinate-start_value)/width)) .* ...
        0.5 .* (1 - tanh((coordinate-end_value)/width));
end
