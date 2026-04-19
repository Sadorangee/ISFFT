clear all; close all; clc;

data = importdata('Poisson_3D_IO.dat');

Lx = 4; Ly = 4; Lz = 4;
nx = 128; ny = 128; nz = 128;

dx = Lx/nx; dy = Ly/ny; dz = Lz/nz;

x = 0:dx:Lx-dx;  y = 0:dy:Ly-dy;  z = 0:dz:Lz-dz;

% Exact solution:
% u(x,y,z) = sin(2*pi*z/Lz) * ( sin(2*pi*x/Lx) + cos(4*pi*y/Ly) )
%
% Periodic Poisson problem:
% Delta u = f
%
% with
% f(x,y,z) = -sin(2*pi*z/Lz) * [ ((2*pi/Lx)^2 + (2*pi/Lz)^2) * sin(2*pi*x/Lx)
%                              + ((4*pi/Ly)^2 + (2*pi/Lz)^2) * cos(4*pi*y/Ly) ]

sx = sin(2*pi*x(:)/Lx);      % [nx,1]
cy = cos(4*pi*y(:)/Ly);      % [ny,1]
sz = sin(2*pi*z(:)/Lz);      % [nz,1]

A = sx * sz.';               % [nx,nz] -> sin(2pi x/Lx)*sin(2pi z/Lz)
B = cy * sz.';               % [ny,nz] -> cos(4pi y/Ly)*sin(2pi z/Lz)
u = permute(A,[1 3 2]) + permute(B,[3 1 2]);   % [nx,ny,nz]

pp = reshape(data(:), [nx,ny,nz]);

figure
plot(z, squeeze(u(nx/2,ny/2,:)), 'K-','LineWidth',2); hold on
plot(z, squeeze(pp(nx/2,ny/2,:)),'o','LineWidth',2,'Color','#C2817C');
xlabel('$z$','Interpreter','latex','FontSize',16); ylabel('$\phi(x_0,y_0,z)$','Interpreter','latex','FontSize',16);
legend('Theory','$SPLASH$','Location',      'northeast', ...
    'Orientation',   'vertical', ...
    'Box',           'off', ...
    'FontSize',      11, ...
    'FontName',      'Times New Roman','Interpreter','latex')

ax = (2*pi/Lx)^2 + (2*pi/Lz)^2;
ay = (4*pi/Ly)^2 + (2*pi/Lz)^2;

AX = ax * (sx * sz.');       % [nx,nz]
BY = ay * (cy * sz.');       % [ny,nz]
f  = -( permute(AX,[1 3 2]) + permute(BY,[3 1 2]) );   % [nx,ny,nz]

kx = 2*pi/Lx*[0:nx/2-1, -nx/2:-1]';   % [nx,1]
ky = 2*pi/Ly*[0:ny/2-1, -ny/2:-1]';   % [ny,1]
kz = 2*pi/Lz*[0:nz/2-1, -nz/2:-1]';   % [nz,1]

Phi = fftn(pp);

LapPhiHat = - bsxfun(@times, reshape(kx.^2,[nx,1,1]), Phi);
LapPhiHat = LapPhiHat - bsxfun(@times, reshape(ky.^2,[1,ny,1]), Phi);
LapPhiHat = LapPhiHat - bsxfun(@times, reshape(kz.^2,[1,1,nz]), Phi);

Lpp_spec = ifftn(LapPhiHat, 'symmetric');

R = Lpp_spec - f;

dV   = dx*dy*dz;
resL2 = sqrt( sum(R(:).^2) * dV );
fprintf('||Δφ - f||_L2 = %.6e\n', resL2);

pp0 = pp - mean(pp(:));
MAE = mean( abs(pp0(:) - u(:)) );
fprintf('MAE(solution vs exact) = %.6e\n', MAE);

set(gcf, 'Position', [100, 100, 500, 400]);
set(gca, 'LineWidth', 1.3);

% print('profile','-dpng','-r300');
