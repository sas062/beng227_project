% Perplexity code...

clc; clear; close all;

load('islet_results.mat');   % expects t, y, p

% -----------------------------
% Indices
% -----------------------------
ia = 1:p.Na;
ib = p.Na + (1:p.Nb);
iG = p.Na + p.Nb + (1:p.Nx*p.Ny);

theta_a = y(:, ia);
theta_b = y(:, ib);
G_all = y(:, iG);

% -----------------------------
% User settings
% -----------------------------
window_fraction = 0.2;     % analyze last 20% of simulation
tol_rel_G = 1e-3;          % relative change tolerance for G
tol_rel_freq = 1e-3;       % relative change tolerance for frequencies
tol_std_freq = 1e-3;       % tolerance for variation within final window

% -----------------------------
% Final analysis window
% -----------------------------
Nt = length(t);
k0 = max(2, floor((1 - window_fraction)*Nt));
tw = t(k0:end);

theta_a_w = theta_a(k0:end, :);
theta_b_w = theta_b(k0:end, :);
G_w = G_all(k0:end, :);

% -----------------------------
% 1) Glucagon field settling
% -----------------------------
G_start = G_w(1, :);
G_end   = G_w(end, :);

rel_change_G = norm(G_end - G_start) / max(norm(G_end), 1e-12);

% frame-to-frame relative changes in final window
rel_step_G = zeros(length(tw)-1, 1);
for k = 1:length(tw)-1
    num = norm(G_w(k+1,:) - G_w(k,:));
    den = max(norm(G_w(k+1,:)), 1e-12);
    rel_step_G(k) = num / den;
end
mean_rel_step_G = mean(rel_step_G);
max_rel_step_G  = max(rel_step_G);

% -----------------------------
% 2) Instantaneous angular velocities
% -----------------------------
theta_a_u = unwrap(theta_a_w);
theta_b_u = unwrap(theta_b_w);

dtw = diff(tw);

omega_a_inst = diff(theta_a_u, 1, 1) ./ dtw;
omega_b_inst = diff(theta_b_u, 1, 1) ./ dtw;

mean_omega_a_t = mean(omega_a_inst, 2);
mean_omega_b_t = mean(omega_b_inst, 2);

std_omega_a_t = std(omega_a_inst, 0, 2);
std_omega_b_t = std(omega_b_inst, 0, 2);

rel_change_mean_omega_a = abs(mean_omega_a_t(end) - mean_omega_a_t(1)) / max(abs(mean_omega_a_t(end)), 1e-12);
rel_change_mean_omega_b = abs(mean_omega_b_t(end) - mean_omega_b_t(1)) / max(abs(mean_omega_b_t(end)), 1e-12);

std_mean_omega_a = std(mean_omega_a_t);
std_mean_omega_b = std(mean_omega_b_t);

% -----------------------------
% 3) Kuramoto order parameters
% -----------------------------
Ra = abs(mean(exp(1i * theta_a_w), 2));
Rb = abs(mean(exp(1i * theta_b_w), 2));

rel_change_Ra = abs(Ra(end) - Ra(1)) / max(abs(Ra(end)), 1e-12);
rel_change_Rb = abs(Rb(end) - Rb(1)) / max(abs(Rb(end)), 1e-12);

std_Ra = std(Ra);
std_Rb = std(Rb);

% -----------------------------
% 4) Fixed-point vs oscillatory regime test
% -----------------------------
fixed_like_G = (rel_change_G < tol_rel_G) && (mean_rel_step_G < tol_rel_G);
freq_settled = ...
    (rel_change_mean_omega_a < tol_rel_freq) && ...
    (rel_change_mean_omega_b < tol_rel_freq) && ...
    (std_mean_omega_a < tol_std_freq) && ...
    (std_mean_omega_b < tol_std_freq);

if fixed_like_G && freq_settled
    if abs(mean_omega_a_t(end)) < tol_std_freq && abs(mean_omega_b_t(end)) < tol_std_freq
        regime = 'Approximate fixed steady state';
    else
        regime = 'Settled oscillatory / phase-locked regime';
    end
else
    regime = 'No clear steady regime yet';
end

% -----------------------------
% 5) Report
% -----------------------------
fprintf('\n===== STEADY-STATE TEST =====\n');
fprintf('Analyzed final %.1f%% of simulation: t in [%.2f, %.2f] s\n', ...
    100*window_fraction, tw(1), tw(end));

fprintf('\nGlucagon field metrics:\n');
fprintf('  Relative change over final window      = %.3e\n', rel_change_G);
fprintf('  Mean relative step change in final win = %.3e\n', mean_rel_step_G);
fprintf('  Max  relative step change in final win = %.3e\n', max_rel_step_G);

fprintf('\nAlpha frequency metrics:\n');
fprintf('  Mean omega_a at end                    = %.6e\n', mean_omega_a_t(end));
fprintf('  Relative drift in mean omega_a         = %.3e\n', rel_change_mean_omega_a);
fprintf('  Std of mean omega_a in final window    = %.3e\n', std_mean_omega_a);

fprintf('\nBeta frequency metrics:\n');
fprintf('  Mean omega_b at end                    = %.6e\n', mean_omega_b_t(end));
fprintf('  Relative drift in mean omega_b         = %.3e\n', rel_change_mean_omega_b);
fprintf('  Std of mean omega_b in final window    = %.3e\n', std_mean_omega_b);

fprintf('\nSynchronization metrics:\n');
fprintf('  Final R_a                              = %.4f\n', Ra(end));
fprintf('  Final R_b                              = %.4f\n', Rb(end));
fprintf('  Relative drift in R_a                  = %.3e\n', rel_change_Ra);
fprintf('  Relative drift in R_b                  = %.3e\n', rel_change_Rb);
fprintf('  Std of R_a in final window             = %.3e\n', std_Ra);
fprintf('  Std of R_b in final window             = %.3e\n', std_Rb);

fprintf('\nClassification:\n');
fprintf('  %s\n', regime);

% -----------------------------
% 6) Plots
% -----------------------------
figure('Color','w', 'Position', [100 100 1100 800]);

subplot(2,2,1);
plot(tw(2:end), mean_omega_a_t, 'g', 'LineWidth', 1.5); hold on;
plot(tw(2:end), mean_omega_b_t, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Mean angular velocity');
title('Mean instantaneous frequencies');
legend('Alpha', 'Beta', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(tw, Ra, 'g', 'LineWidth', 1.5); hold on;
plot(tw, Rb, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Order parameter R');
title('Synchronization over final window');
legend('Alpha', 'Beta', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(tw(2:end), rel_step_G, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Relative step change');
title('Glucagon field change per time step');
grid on;

subplot(2,2,4);
plot(tw, mean(G_w, 2), 'b', 'LineWidth', 1.5); hold on;
plot(tw, std(G_w, 0, 2), 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Value');
title('Mean and std of glucagon field');
legend('Mean(G)', 'Std(G)', 'Location', 'best');
grid on;
