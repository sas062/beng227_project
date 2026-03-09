clc; clear; close all;

% Solver parameters
tspan = [0 720];

% Define constants
p.w_a = 2*pi / 40;
p.tau = 0.05;
p.DG = 100;      % placeholder diffusivity
p.sG = 200;       % placeholder secretion scaling factor
p.sigmaG = 10;   % placeholder secretion Gaussian SD. ~3*sigmaG = max secretion radius
p.lI = 20;       % placeholder interaction length (for insulin secretion to alpha cells)

% Coupling Coefficients
p.K_ab = 0.52 * (0.2*pi);
p.K_ba = -0.66 * (0.2*pi);
p.K_bb = 0.08;

% Original Secretion Response Functions
% Ren et al Equations 1-3
p.fs = @(th) (1-cos(th)) / 2;
p.fra = @(th) sin(th);
p.frb = @(th) 1;

N = 100; % TOTAL NUMBER OF CELLS
Na = floor(N*0.15); % NUMBER OF ALPHA CELLS (based off ~ 15% value from literature)
Nb = N-Na; % NUMBER OF BETA CELLS
p.Nb = Nb; % NUMBER OF BETA CELLS
p.Na = Na; % NUMBER OF ALPHA CELLS

% Heterogeneity of Beta-cell intrinsic frequency
% p.w_b : column vector of heterogenous beta frequencies
p.wb_mean = 2*pi/360; % mean frequency (from Ren)
p.sigma_b = 0.5; % Standard deviation for heterogeneity
rng(1); % reproducible rand seed

% Parameters from distribution fitting (gamma function)
a = 2.1024;
b = 1.98715;
p.w_b = 2*pi./(random('gamma',a,b,[p.Nb 1])*60);
%p.w_b = p.wb_mean * (1 + p.sigma_b * randn(p.Nb,1));


% Initial alpha and beta phases
% Could also just set to zeros
theta_a0 = 2*pi * rand(p.Na,1);
theta_b0 = 2*pi * rand(p.Nb,1);

% 2D spatial grid for glucagon
p.Lx = 400; % in μm
p.Ly = 400;
p.Nx = 81; % approx 5 μm resolution (40 intervals)
p.Ny = 81;

p.x = linspace(0, p.Lx, p.Nx);
p.y = linspace(0, p.Ly, p.Ny);
[p.X, p.Y] = meshgrid(p.x, p.y);

p.dx = p.x(2) - p.x(1);
p.dy = p.y(2) - p.y(1);

% Create Islet Structure
Ri = 10;   % um
thr = 100; % number of errors before increasing islet radius
r = 5;     % radius of cell

[x,y,Rf] = circle_packing(Ri,r,thr,N);
arrangement = "toward_outer";
[xa,ya,xb,yb] = assign_cell_type(x,y,Na,Nb,arrangement);
xa = xa';
xb = xb';
ya = ya';
yb = yb';

% Circular islet fit inside grid
p.cx = p.Lx / 2;
p.cy = p.Ly / 2;
p.R_islet = Rf;
p.xa = p.cx + xa;
p.ya = p.cy + ya;
p.xb = p.cx + xb;
p.yb = p.cy + yb;

thr_r = 20;
NN = get_NN(xb,yb,thr_r);
p.NN = NN;

plot_structure(p) % Visualize Islet Structure
% % Circular islet fit inside grid
% p.cx = p.Lx / 2;
% p.cy = p.Ly / 2;
% p.R_islet = 0.45 * min(p.Lx,p.Ly);
% 
% rb = p.R_islet * sqrt(0.8 * rand(p.Nb,1));
% phi_b = 2*pi * rand(p.Nb,1);
% p.xb = p.cx + rb .* cos(phi_b);
% p.yb = p.cy + rb .* sin(phi_b);
% 
% ra = p.R_islet * (0.75 + 0.25 * rand(p.Na,1));
% phi_a = 2*pi * rand(p.Na,1);
% p.xa = p.cx + ra .* cos(phi_a);
% p.ya = p.cy + ra .* sin(phi_a);

% Initial Conditions
p.diffusion = 1;
switch p.diffusion
    case 1
        G0 = zeros(p.Ny, p.Nx);
        y0 = [theta_a0; theta_b0; G0(:)];
        tic
        [t, y] = ode15s(@(t,y) rhs_islet(t,y,p), tspan, y0);
        save('islet_results.mat', 't', 'y', 'p', 'theta_a0', 'theta_b0', 'G0');
        toc
    case 0 % Runs model without diffusion, only beta-beta coupling -> need this to get an idea of what are K_bb should be
        G0 = zeros(p.Ny, p.Nx);
        y0 = [theta_a0; theta_b0];
        tic
        [t, y] = ode15s(@(t,y) rhs_islet_no_glucagon(t,y,p), tspan, y0);
        toc
        save('islet_results.mat', 't', 'y', 'p', 'theta_a0', 'theta_b0','G0');      
end




%% Just for quick looks
% Extract phase histories
theta_a = y(:, 1:p.Na);
theta_b = y(:, p.Na + (1:p.Nb));
theta_a = mod(theta_a,2*pi);
theta_b = mod(theta_b,2*pi);

figure(3);

%subplot(1,2,1);
plot(t, theta_a, 'LineWidth', 1);
xlabel('Time');
ylabel('\theta_\alpha (rad)');
title('Alpha-cell phases');
grid on;


figure(4);
%subplot(1,2,2);
plot(t, theta_b, 'LineWidth', 1);
xlabel('Time');
ylabel('\theta_\beta (rad)');
title('Beta-cell phases');
grid on;

function plot_structure(p)
r = 5;
xa = p.xa;
xb = p.xb;
ya = p.ya;
yb = p.yb;
Na = p.Na;
Nb = p.Nb;
NN = p.NN;
Rf = p.R_islet;
figure;
for i=1:Na
    rectangle('Position',[xa(i)-r ya(i)-r 2*r 2*r],'Curvature',[1,1],'FaceColor','g','EdgeColor','g');
end

for i=1:Nb
    rectangle('Position',[xb(i)-r yb(i)-r 2*r 2*r],'Curvature',[1,1],'FaceColor', 'r', 'EdgeColor','r');
end
rectangle('Position',[-Rf+p.cx -Rf+p.cy 2*Rf 2*Rf],'Curvature',[1,1],'EdgeColor','k')

for i = 1:Nb
    xi = xb(i);
    yi = yb(i);
    NNi = NN(i,:);
    [~,cells,~] = find(NNi);
    for j = 1:numel(cells)
        xj = xb(cells(j));
        yj = yb(cells(j));
        line([xi xj], [yi yj],'Color', 'black','LineWidth',1)
    end    
end

end