function p = ren2022_default_params()
%REN2022_DEFAULT_PARAMS Parameter set for Ren et al. (Nat Commun, 2022).
%   This model follows the coupled phase-oscillator formulation used in the
%   paper to describe global phase locking between pancreatic alpha and beta
%   cells through paracrine signaling.
%
%   Equations:
%     d theta_alpha / dt = omega_alpha + K_beta_to_alpha * f_s(theta_beta) * f_r_alpha(theta_alpha)
%     d theta_beta  / dt = omega_beta  + K_alpha_to_beta * f_s(theta_alpha) * f_r_beta(theta_beta)
%
%   where f_s is a secretion-like activity profile and f_r_* are
%   phase-response profiles.

% Intrinsic cycle periods (seconds). Adjust if needed for your data.
p.T_alpha = 50;   % alpha intrinsic period
p.T_beta  = 240;  % beta intrinsic period

% Intrinsic angular frequencies (rad/s)
p.omega_alpha = 2*pi / p.T_alpha;
p.omega_beta  = 2*pi / p.T_beta;

% Coupling gains.
% K_alpha_to_beta > 0 : glucagon-like stimulation of beta cells
% K_beta_to_alpha < 0 : insulin/somatostatin-like inhibition of alpha cells
p.K_alpha_to_beta = 0.25;
p.K_beta_to_alpha = -0.45;

% Secretion pulse shape parameters
p.secretion_power = 3;   % sharpness of secretion pulse

% Phase-response strength parameters
p.response_alpha_amp = 1.0;
p.response_beta_amp  = 1.0;

% Simulation settings
p.tspan = [0, 3600];   % 1 hour
p.theta0 = [0.3; 3.5]; % initial phases [alpha; beta]
end
