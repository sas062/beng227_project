function dtheta = ren2022_rhs(~, theta, p)
%REN2022_RHS Right-hand side for the Ren et al. 2022 phase model.
%   theta(1) = alpha-cell phase
%   theta(2) = beta-cell phase

theta_a = theta(1);
theta_b = theta(2);

% Secretion-like profile: active mainly in the rising/active phase.
fs_a = secretion_profile(theta_a, p.secretion_power);
fs_b = secretion_profile(theta_b, p.secretion_power);

% Phase-response profiles (sign indicates stimulatory or inhibitory effect).
fr_a = p.response_alpha_amp * (1 - cos(theta_a));
fr_b = p.response_beta_amp  * (1 - cos(theta_b));

% Coupled phase equations.
dtheta_a = p.omega_alpha + p.K_beta_to_alpha * fs_b * fr_a;
dtheta_b = p.omega_beta  + p.K_alpha_to_beta * fs_a * fr_b;

dtheta = [dtheta_a; dtheta_b];
end

function y = secretion_profile(theta, q)
%SECRETION_PROFILE Positive, periodic secretion pulse from phase.
%   y in [0,1], peaking near theta = pi/2 and zero in silent half-cycle.
y = max(0, sin(theta));
y = y.^q;
end
