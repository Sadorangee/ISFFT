%% Settings. Keep these consistent with Poisson_solver_UI.in.
boundary = [0 0 0];       % 0 periodic, 1 homogeneous Dirichlet, 2 homogeneous Neumann
dimension = 3;            % Active directions are x, then y, then z.
grid_size = [128 128 128];
domain_size = [4 4 4];

%% Read the solution and construct the directional modes.
validateattributes(dimension, {'numeric'}, {'scalar','integer','>=',1,'<=',3});
validateattributes(boundary, {'numeric'}, {'vector','numel',3,'integer','finite'});
validateattributes(grid_size, {'numeric'}, {'vector','numel',3,'integer','positive','finite'});
validateattributes(domain_size, {'numeric'}, {'vector','numel',3,'positive','finite'});
boundary = boundary(:).';
grid_size = grid_size(:).';
domain_size = domain_size(:).';
assert(all(ismember(boundary(1:dimension), 0:2)), 'Active boundary codes must be 0, 1, or 2.');
assert(all(grid_size(dimension+1:3) == 1), 'Inactive grid sizes must be one.');

data_directory = '.';
data_file = fullfile(data_directory, sprintf('Poisson_%dD_IO.dat', dimension));
assert(isfile(data_file), 'Solution file does not exist: %s', data_file);
data = load(data_file, '-ascii');
assert(numel(data) == prod(grid_size), 'Data count differs from the configured grid size.');
assert(all(isfinite(data(:))), 'Solution data contain nonfinite values.');
numerical_solution = reshape(data, grid_size);

coordinate = {0, 0, 0};
mode_value = {1, 1, 1};
quadrature_weight = {1, 1, 1};
wave_number = zeros(1, 3);

for axis = 1:dimension
    node_count = grid_size(axis);
    domain_length = domain_size(axis);
    wave_number(axis) = pi/domain_length;
    if boundary(axis) == 0
        grid_spacing = domain_length/node_count;
        wave_number(axis) = 2*pi/domain_length;
        quadrature_weight{axis} = grid_spacing*ones(node_count, 1);
    else
        assert(node_count >= 2, 'Nonperiodic axes require at least two grid points.');
        grid_spacing = domain_length/(node_count-1);
        quadrature_weight{axis} = grid_spacing*ones(node_count, 1);
        quadrature_weight{axis}([1 end]) = grid_spacing/2;
    end
    coordinate{axis} = (0:node_count-1)'*grid_spacing;
    if boundary(axis) == 1
        mode_value{axis} = sin(wave_number(axis)*coordinate{axis});
    else
        mode_value{axis} = cos(wave_number(axis)*coordinate{axis});
    end
end

% The same continuous analytical solution is used for all discretization schemes.
analytical_solution = reshape(mode_value{1}, [], 1, 1) ...
    .* reshape(mode_value{2}, 1, [], 1) .* reshape(mode_value{3}, 1, 1, []);
volume_weight = reshape(quadrature_weight{1}, [], 1, 1) ...
    .* reshape(quadrature_weight{2}, 1, [], 1) .* reshape(quadrature_weight{3}, 1, 1, []);

% P/N-only problems determine the solution up to an additive constant.
mean_adjustment = 0;
if all(boundary(1:dimension) ~= 1)
    mean_adjustment = sum(volume_weight(:).*(numerical_solution(:)-analytical_solution(:))) ...
        /sum(volume_weight(:));
    numerical_solution = numerical_solution-mean_adjustment;
end
solution_error = numerical_solution-analytical_solution;
error_Linf = max(abs(solution_error(:)));
error_L2 = sqrt(sum(volume_weight(:).*solution_error(:).^2));

%% Evaluate the residual with the discretization used by the solver.
input_text = fileread(fullfile(data_directory, 'Poisson_solver_UI.in'));
scheme_entry = regexp(input_text, ...
    '(?m)^\s*#\s*Scheme_type[^\r\n]*\r?\n([^\r\n]*)', 'tokens', 'once');
assert(~isempty(scheme_entry), 'Cannot find Scheme_type in Poisson_solver_UI.in.');
discretization = sscanf(scheme_entry{1}, '%f');
validateattributes(discretization, {'numeric'}, {'scalar','integer','>=',0,'<=',6});

laplacian = zeros(grid_size);
for axis = 1:dimension
    laplacian = laplacian + axis_laplacian(numerical_solution, axis, ...
        boundary(axis), domain_size(axis), discretization);
end
source_term = -sum(wave_number.^2)*analytical_solution;
residual_field = laplacian-source_term;

% Dirichlet endpoints carry the boundary equation phi = 0.
for axis = 1:dimension
    if boundary(axis) == 1
        boundary_index = {':', ':', ':'};
        boundary_index{axis} = [1 grid_size(axis)];
        residual_field(boundary_index{:}) = numerical_solution(boundary_index{:});
    end
end
residual_Linf = max(abs(residual_field(:)));

fprintf('Linf      %.16e\n', error_Linf);
fprintf('L2        %.16e\n', error_L2);
fprintf('Residuals %.16e\n', residual_Linf);

%% Plot one profile. Other directions are fixed near a mode maximum.
profile_axis = 1;
axis_name = {'x', 'y', 'z'};
profile_index = {1, 1, 1};
profile_factor = 1;
for axis = 1:dimension
    if axis == profile_axis
        continue
    end
    [~, node] = max(abs(mode_value{axis}));
    profile_index{axis} = node;
    profile_factor = profile_factor*mode_value{axis}(node);
end
profile_index{profile_axis} = ':';
numerical_profile = numerical_solution(profile_index{:});
profile_coordinate = linspace(0, domain_size(profile_axis), 1001);
if boundary(profile_axis) == 1
    analytical_profile = profile_factor*sin(wave_number(profile_axis)*profile_coordinate);
else
    analytical_profile = profile_factor*cos(wave_number(profile_axis)*profile_coordinate);
end

figure('Color', 'w', 'Renderer', 'painters');
plot(profile_coordinate, analytical_profile, 'k-', 'LineWidth', 1.5);
hold on
plot(coordinate{profile_axis}, numerical_profile(:), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 5, 'Color', [0.1 0.4 0.8]);
hold off
grid on
xlabel(axis_name{profile_axis});
ylabel('\phi');
title(sprintf('%dD Poisson solution, %s profile', dimension, axis_name{profile_axis}));
legend('Analytical solution', 'Numerical solution', 'Location', 'best');

function result = axis_laplacian(field, axis, boundary_type, domain_length, discretization)
    node_count = size(field, axis);
    index = {':', ':', ':'};
    if boundary_type == 0
        grid_spacing = domain_length/node_count;
        wave_scale = 2*pi/domain_length;
    else
        grid_spacing = domain_length/(node_count-1);
        wave_scale = pi/domain_length;
        if boundary_type == 1
            index{axis} = [1 node_count];
            field(index{:}) = 0;
        end
        index{axis} = node_count-1:-1:2;
        reflected = field(index{:});
        if boundary_type == 1
            reflected = -reflected;
        end
        % Odd extension gives the sine basis, even extension the cosine basis.
        field = cat(axis, field, reflected);
    end

    transform_size = size(field, axis);
    mode = (0:transform_size-1)';
    mode = min(mode, transform_size-mode);
    wave = wave_scale*mode;
    if discretization == 0
        lambda = wave.^2;
    else
        % Columns are alpha, beta, a, b, c from ISFFT_Discrete_Operators.f90.
        coefficients = [ ...
            0, 0, 1, 0, 0; ...
            0, 0, 4/3, -1/3, 0; ...
            1/10, 0, 6/5, 0, 0; ...
            0, 0, 3/2, -3/5, 1/10; ...
            2/11, 0, 12/11, 3/11, 0; ...
            0.502092570454933, 0.055691658873965, ...
            0.215649631004306, 1.723321721695239, 0.176597105958252];
        alpha = coefficients(discretization,1);
        beta = coefficients(discretization,2);
        a = coefficients(discretization,3);
        b = coefficients(discretization,4);
        c = coefficients(discretization,5);
        theta = wave*grid_spacing;
        lambda = (4*a*sin(theta/2).^2 + b*sin(theta).^2 + (4*c/9)*sin(3*theta/2).^2) ...
            ./ (grid_spacing^2*(1 + 2*alpha*cos(theta) + 2*beta*cos(2*theta)));
    end

    scale_size = [1 1 1];
    scale_size(axis) = transform_size;
    spectrum = fft(field, [], axis);
    result = real(ifft(-reshape(lambda, scale_size).*spectrum, [], axis));
    index{axis} = 1:node_count;
    result = result(index{:});
end
