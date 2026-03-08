clc; clear; close all;

% Solver parameters
tspan = [0 3600];

% Define constants
p.w_a = 2*pi / 40;
p.tau = 0.5;
p.DG = 100;      % placeholder diffusivity
p.sG = 10;       % placeholder secretion scaling factor
p.sigmaG = 25;   % placeholder secretion Gaussian SD. ~3*sigmaG = max secretion radius
p.lI = 20;       % placeholder interaction length (for insulin secretion to alpha cells)

% Coupling Coefficients
p.K_ab = 0.52 * (0.2*pi);
p.K_ba = -0.66 * (0.2*pi);

% Original Secretion Response Functions
% Ren et al Equations 1-3
p.fs = @(th) (1-cos(th)) / 2;
p.fra = @(th) sin(th);
p.frb = @(th) 1;

% Heterogeneity of Beta-cell intrinsic frequency
% p.w_b : column vector of heterogenous beta frequencies
p.Nb = 40; % NUMBER OF BETA CELLS
p.wb_mean = 2*pi/360; % mean frequency (from Ren)
p.sigma_b = 0.1; % Standard deviation for heterogeneity
rng(1); % reproducible rand seed
p.w_b = p.wb_mean * (1 + p.sigma_b * randn(p.Nb,1));

p.Na = 20; % NUMBER OF ALPHA CELLS

% Initial alpha and beta phases
% Could also just set to zeros
theta_a0 = 2*pi * rand(p.Na,1);
theta_b0 = 2*pi * rand(p.Nb,1);

% 2D spatial grid for glucagon
p.Lx = 200; % in μm
p.Ly = 200;
p.Nx = 41; % approx 5 μm resolution (40 intervals)
p.Ny = 41;

p.x = linspace(0, p.Lx, p.Nx);
p.y = linspace(0, p.Ly, p.Ny);
[p.X, p.Y] = meshgrid(p.x, p.y);

p.dx = p.x(2) - p.x(1);
p.dy = p.y(2) - p.y(1);

% 2D cell/grid layout with Beta cells interior, Alpha cells peripheral

% Circular islet fit inside grid
p.cx = p.Lx / 2;
p.cy = p.Ly / 2;
p.R_islet = 0.45 * min(p.Lx,p.Ly);

rb = p.R_islet * sqrt(0.8 * rand(p.Nb,1));
phi_b = 2*pi * rand(p.Nb,1);
p.xb = p.cx + rb .* cos(phi_b);
p.yb = p.cy + rb .* sin(phi_b);

ra = p.R_islet * (0.75 + 0.25 * rand(p.Na,1));
phi_a = 2*pi * rand(p.Na,1);
p.xa = p.cx + ra .* cos(phi_a);
p.ya = p.cy + ra .* sin(phi_a);

% Initial Conditions
G0 = zeros(p.Ny, p.Nx);
y0 = [theta_a0; theta_b0; G0(:)];

% SOLVE
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t, y] = ode15s(@(t,y) rhs_islet(t,y,p), tspan, y0);

save('islet_results.mat', 't', 'y', 'p', 'theta_a0', 'theta_b0', 'G0');