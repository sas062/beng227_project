function out = run_ren2022_demo(userParams)
%RUN_REN2022_DEMO Simulate alpha-beta phase locking model from Ren et al.
%   OUT = RUN_REN2022_DEMO() runs with default parameters.
%   OUT = RUN_REN2022_DEMO(P) overrides fields in the default parameter
%   struct with matching fields from P.
%
%   Plot rationale:
%   1) Wrapped phases: verify oscillatory activity and timing relation.
%   2) Relative phase over time: direct signature of phase locking
%      (plateau/confinement indicates locking).
%   3) Relative phase histogram: reveals preferred phase difference.
%   4) Cycle counts + winding ratio: checks global locking ratio n_alpha/n_beta.

p = ren2022_default_params();
if nargin > 0 && ~isempty(userParams)
    names = fieldnames(userParams);
    for k = 1:numel(names)
        p.(names{k}) = userParams.(names{k});
    end
end

odefun = @(t, y) ren2022_rhs(t, y, p);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
[t, theta] = ode45(odefun, p.tspan, p.theta0, opts);

theta_a = theta(:,1);
theta_b = theta(:,2);

% Wrapped phases for visualization
phi_a = mod(theta_a, 2*pi);
phi_b = mod(theta_b, 2*pi);

% Relative phase wrapped to [-pi, pi]
rel = angle(exp(1i*(phi_a - phi_b)));

% Estimate locking ratio via winding numbers.
na = (theta_a(end) - theta_a(1)) / (2*pi);
nb = (theta_b(end) - theta_b(1)) / (2*pi);
ratio = na / nb;

% Instantaneous angular velocities (numerical derivatives)
omega_a_inst = gradient(theta_a, t);
omega_b_inst = gradient(theta_b, t);
omega_ratio_inst = omega_a_inst ./ max(omega_b_inst, eps);

out = struct();
out.t = t;
out.theta_alpha = theta_a;
out.theta_beta = theta_b;
out.phi_alpha = phi_a;
out.phi_beta = phi_b;
out.relative_phase = rel;
out.winding_alpha = na;
out.winding_beta = nb;
out.winding_ratio = ratio;
out.omega_alpha_inst = omega_a_inst;
out.omega_beta_inst = omega_b_inst;
out.omega_ratio_inst = omega_ratio_inst;
out.params = p;

fprintf('Estimated winding ratio n_alpha/n_beta = %.3f\n', ratio);

% Plots
figure('Name', 'Ren et al. 2022 phase-oscillator model', 'Color', 'w');

subplot(2,2,1);
plot(t, phi_a, 'r', 'LineWidth', 1.1); hold on;
plot(t, phi_b, 'b', 'LineWidth', 1.1);
ylabel('wrapped phase (rad)');
legend({'\alpha', '\beta'}, 'Location', 'best');
title('Wrapped phases');
grid on;

subplot(2,2,2);
plot(t, rel, 'k', 'LineWidth', 1.1);
ylabel('\phi_\alpha - \phi_\beta (rad)');
title('Relative phase dynamics');
grid on;

subplot(2,2,3);
nbins = 36;
histogram(rel, nbins, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
xlabel('\phi_\alpha - \phi_\beta (rad)');
ylabel('count');
title('Relative phase distribution');
grid on;

subplot(2,2,4);
plot(t, theta_a/(2*pi), 'r', 'LineWidth', 1.0); hold on;
plot(t, theta_b/(2*pi), 'b', 'LineWidth', 1.0);
ylabel('cycle count');
xlabel('time (s)');
title(sprintf('Cycle winding ratio n_\\alpha/n_\\beta = %.3f', ratio));
grid on;

end
