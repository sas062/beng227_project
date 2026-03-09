% Perplexity code...

% Plots

clc; clear; close all;
load('islet_results.mat');

% Initial and final fields
S0 = glucagon_secretion(theta_a0, p);
G_final = reshape(y(end, p.Na + p.Nb + 1:end), p.Ny, p.Nx);

% Figure 1: surf-style comparison
figure(1);

subplot(1,2,1);
surf(p.X, p.Y, S0);
shading interp;
view(2);
axis equal tight;
cb1 = colorbar;
cb1.Color = 'k';
hold on;
z1 = max(S0(:)) + 0.05*max(S0(:));
hA1 = scatter3(p.xa, p.ya, z1*ones(size(p.xa)), 60, 'g', 'filled', 'MarkerEdgeColor', 'k');
hB1 = scatter3(p.xb, p.yb, z1*ones(size(p.xb)), 60, 'r', 'filled', 'MarkerEdgeColor', 'k');
title('Initial Glucagon Secretion');
xlabel('x (\mum)');
ylabel('y (\mum)');

subplot(1,2,2);
surf(p.X, p.Y, G_final);
shading interp;
view(2);
axis equal tight;
cb2 = colorbar;
cb2.Color = 'k';
hold on;
z1 = max(S0(:)) + 0.05*max(S0(:));
hA2 = scatter3(p.xa, p.ya, z1*ones(size(p.xa)), 60, 'g', 'filled', 'MarkerEdgeColor', 'k');
hB2 = scatter3(p.xb, p.yb, z1*ones(size(p.xb)), 60, 'r', 'filled', 'MarkerEdgeColor', 'k');
title(sprintf('Final Glucagon Concentration, t = %.2f', t(end)));
xlabel('x (\mum)');
ylabel('y (\mum)');
legend([hA2, hB2], {'Alpha cells', 'Beta cells'}, 'Location', 'best');

% Figure 2: imagesc comparison
figure(2);

subplot(1,2,1);
imagesc(p.x, p.y, S0);
axis xy equal tight;
cb3 = colorbar;
cb3.Color = 'k';
hold on;
hA3 = plot(p.xa, p.ya, 'wo', 'MarkerFaceColor', 'g');
hB3 = plot(p.xb, p.yb, 'ko', 'MarkerFaceColor', 'r');
title('Initial Glucagon Secretion');
xlabel('x (\mum)');
ylabel('y (\mum)');

subplot(1,2,2);
imagesc(p.x, p.y, G_final);
axis xy equal tight;
cb4 = colorbar;
cb4.Color = 'k';
hold on;
hA4 = plot(p.xa, p.ya, 'wo', 'MarkerFaceColor', 'g');
hB4 = plot(p.xb, p.yb, 'ko', 'MarkerFaceColor', 'r');
title(sprintf('Final Glucagon Concentration, t = %.2f', t(end)));
xlabel('x (\mum)');
ylabel('y (\mum)');
legend([hA4, hB4], {'Alpha cells', 'Beta cells'}, 'Location', 'best');

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

figure(3);

subplot(1,2,1);
plot(t, theta_a, 'LineWidth', 1);
xlabel('Time');
ylabel('\theta_\alpha (rad)');
title('Alpha-cell phases');
grid on;

subplot(1,2,2);
plot(t, theta_b, 'LineWidth', 1);
xlabel('Time');
ylabel('\theta_\beta (rad)');
title('Beta-cell phases');
grid on;

theta_a_unwrapped = unwrap(theta_a);
theta_b_unwrapped = unwrap(theta_b);

figure(4);

subplot(1,2,1);
plot(t, theta_a_unwrapped, 'LineWidth', 1);
xlabel('Time');
ylabel('Unwrapped \theta_\alpha (rad)');
title('Alpha-cell phases (unwrapped)');
grid on;

subplot(1,2,2);
plot(t, theta_b_unwrapped, 'LineWidth', 1);
xlabel('Time');
ylabel('Unwrapped \theta_\beta (rad)');
title('Beta-cell phases (unwrapped)');
grid on;

mean_theta_a = mean(theta_a_unwrapped, 2);
mean_theta_b = mean(theta_b_unwrapped, 2);

figure(5);
plot(t, mean_theta_a, 'g', 'LineWidth', 2); hold on;
plot(t, mean_theta_b, 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('Mean unwrapped phase (rad)');
title('Population mean phases');
legend('Alpha mean phase', 'Beta mean phase', 'Location', 'best');
grid on;

figure(6);
plot(t, mean_theta_a - mean_theta_b, 'b', 'LineWidth', 2);
xlabel('Time');
ylabel('\Delta phase (rad)');
title('Mean alpha-beta phase difference');
grid on;

% Choose representative cells
ia_plot = 1;
ib_plot = 1;

theta_a = y(:, ia_plot);
theta_b = y(:, p.Na + ib_plot);

figure(7);
plot(t, mod(theta_a, 2*pi), 'g', 'LineWidth', 2); hold on;
plot(t, mod(theta_b, 2*pi), 'r', 'LineWidth', 2);

xlabel('Time');
ylabel('\theta');
title('Hybrid model phase comparison');
legend('\theta_a', '\theta_b', 'Location', 'northeast');
grid on;
xlim([t(1) t(end)]);
ylim([0 2*pi]);


% Better version of fig 7 ?

% Extract all phase histories
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));

% Unwrap to remove artificial 2*pi jumps
theta_a_unwrap = unwrap(theta_a_all);
theta_b_unwrap = unwrap(theta_b_all);

% Compute population-mean trajectories
theta_a_mean = mean(theta_a_unwrap, 2);   % column vector in time
theta_b_mean = mean(theta_b_unwrap, 2);

% Compute distance of each cell trajectory from its population mean
dist_a = sum((theta_a_unwrap - theta_a_mean).^2, 1);
dist_b = sum((theta_b_unwrap - theta_b_mean).^2, 1);

% Pick most representative cells
[~, ia_plot] = min(dist_a);
[~, ib_plot] = min(dist_b);

% Extract those representative trajectories
theta_a_rep = theta_a_all(:, ia_plot);
theta_b_rep = theta_b_all(:, ib_plot);

% Plot in Ren-style wrapped format
figure(8);
plot(t, mod(theta_a_rep, 2*pi), 'g', 'LineWidth', 2); hold on;
plot(t, mod(theta_b_rep, 2*pi), 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('\theta');
title('Representative Hybrid Model Phase Comparison');
legend('\theta_a', '\theta_b', 'Location', 'northeast');
grid on;
xlim([t(1) t(end)]);
ylim([0 2*pi]);


% Anotha one

% Extract all phase histories
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));

% Choose representative alpha cell from trajectory closeness to alpha mean
theta_a_unwrap = unwrap(theta_a_all);
theta_a_mean = mean(theta_a_unwrap, 2);
dist_a = sum((theta_a_unwrap - theta_a_mean).^2, 1);
[~, ia_plot] = min(dist_a);

% Choose beta cell with intrinsic frequency closest to mean beta frequency
[~, ib_plot] = min(abs(p.w_b - mean(p.w_b)));

% Extract selected trajectories
theta_a = theta_a_all(:, ia_plot);
theta_b = theta_b_all(:, ib_plot);

% Plot wrapped phases in Ren-style format
figure(9);
plot(t, mod(theta_a, 2*pi), 'g', 'LineWidth', 2); hold on;
plot(t, mod(theta_b, 2*pi), 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('\theta');
title('Hybrid Model Phase Comparison');
legend('\theta_a', '\theta_b', 'Location', 'northeast');
grid on;
xlim([t(1) t(end)]);
ylim([0 2*pi]);

fprintf('Representative alpha cell: %d\n', ia_plot);
fprintf('Representative beta cell (closest to mean frequency): %d\n', ib_plot);
fprintf('Mean beta frequency: %.6f rad/s\n', mean(p.w_b));
fprintf('Selected beta frequency: %.6f rad/s\n', p.w_b(ib_plot));

%% Do network analysis

load('islet_results.mat');

theta_b = y(:, p.Na + (1:p.Nb));
theta_b = mod(theta_b,2*pi);

st = find(t > 50,1, 'first');
ed = numel(t);
fspan = [st ed];
Kavg_desired = 8;
R = network_analysis(Kavg_desired,theta_b,fspan);