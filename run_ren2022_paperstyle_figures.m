function out = run_ren2022_paperstyle_figures(userParams)
%RUN_REN2022_PAPERSTYLE_FIGURES Paper-style visual checks for Ren et al. model.
%   OUT = RUN_REN2022_PAPERSTYLE_FIGURES() simulates the alpha-beta phase
%   model and draws figure panels arranged to visually match the *style* of
%   the paper's phase-locking summaries (time traces, phase relation, and
%   locking-mode map).
%
%   OUT = RUN_REN2022_PAPERSTYLE_FIGURES(P) overrides default model params.
%
%   This file is intended as a visual validation helper: it uses the same ODE
%   core (`ren2022_rhs`) and derives paper-style diagnostics from the resulting
%   trajectories.

p = ren2022_default_params();
if nargin > 0 && ~isempty(userParams)
    names = fieldnames(userParams);
    for k = 1:numel(names)
        p.(names{k}) = userParams.(names{k});
    end
end

[t, theta] = ode45(@(tt, yy) ren2022_rhs(tt, yy, p), p.tspan, p.theta0, ...
    odeset('RelTol', 1e-7, 'AbsTol', 1e-9));

theta_a = theta(:,1);
theta_b = theta(:,2);
phi_a = mod(theta_a, 2*pi);
phi_b = mod(theta_b, 2*pi);
dphi = angle(exp(1i*(phi_a - phi_b)));

na = (theta_a(end) - theta_a(1)) / (2*pi);
nb = (theta_b(end) - theta_b(1)) / (2*pi);
windingRatio = na / nb;

% Approximate "activation" traces from phase (proxy for Ca2+ activity)
act_a = max(0, sin(phi_a));
act_b = max(0, sin(phi_b));

% Cycle event times (phase wraps) for a raster-like event panel
crossA = find(diff(phi_a) < -pi);
crossB = find(diff(phi_b) < -pi);
event_t_a = t(crossA + 1);
event_t_b = t(crossB + 1);

% --- Figure 1: time-domain and phase-locking summaries ---
figure('Name', 'Ren 2022 paper-style checks: time + phase', 'Color', 'w');

subplot(2,2,1);
plot(t, act_a, 'r', 'LineWidth', 1.2); hold on;
plot(t, act_b, 'b', 'LineWidth', 1.2);
ylabel('activity proxy');
title('A) alpha/beta activity traces (paper-style)');
legend({'\alpha','\beta'}, 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t, dphi, 'k', 'LineWidth', 1.2);
ylabel('\Delta\phi = \phi_\alpha - \phi_\beta (rad)');
title('B) relative phase over time');
grid on;

subplot(2,2,3);
polarhistogram(dphi, 24, 'FaceColor', [0.1 0.1 0.1], 'EdgeColor', 'none');
title('C) preferred locking angle');

subplot(2,2,4);
plot(t, theta_a/(2*pi), 'r', 'LineWidth', 1.0); hold on;
plot(t, theta_b/(2*pi), 'b', 'LineWidth', 1.0);
xlabel('time (s)');
ylabel('cycle count');
title(sprintf('D) winding ratio n_\\alpha/n_\\beta = %.3f', windingRatio));
grid on;

% --- Figure 2: event-raster style and locking-mode map ---
figure('Name', 'Ren 2022 paper-style checks: events + mode map', 'Color', 'w');

subplot(1,2,1);
if ~isempty(event_t_a)
    plot(event_t_a, 1 + 0*event_t_a, 'r|', 'MarkerSize', 10, 'LineWidth', 1.5); hold on;
end
if ~isempty(event_t_b)
    plot(event_t_b, 2 + 0*event_t_b, 'b|', 'MarkerSize', 10, 'LineWidth', 1.5);
end
ylim([0.5 2.5]);
yticks([1 2]);
yticklabels({'\alpha events','\beta events'});
xlabel('time (s)');
title('E) cycle-event raster (phase wraps)');
grid on;

subplot(1,2,2);
[pa, pb, ratioMap] = local_mode_map(p);
imagesc(pb, pa, ratioMap);
set(gca, 'YDir', 'normal');
xlabel('T_\beta (s)');
ylabel('T_\alpha (s)');
title('F) locking-mode map (winding ratio)');
cb = colorbar;
cb.Label.String = 'n_\alpha / n_\beta';
colormap(turbo);
hold on;
plot(p.T_beta, p.T_alpha, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
text(p.T_beta, p.T_alpha, '  current params', 'Color', 'w', 'FontSize', 9);

out = struct();
out.t = t;
out.theta_alpha = theta_a;
out.theta_beta = theta_b;
out.phi_alpha = phi_a;
out.phi_beta = phi_b;
out.relative_phase = dphi;
out.activity_alpha = act_a;
out.activity_beta = act_b;
out.event_t_alpha = event_t_a;
out.event_t_beta = event_t_b;
out.winding_alpha = na;
out.winding_beta = nb;
out.winding_ratio = windingRatio;
out.params = p;

fprintf('Paper-style check complete. Estimated n_alpha/n_beta = %.3f\n', windingRatio);

end

function [TaGrid, TbGrid, ratioMap] = local_mode_map(p)
% Sweep intrinsic periods to create a locking-mode map (Arnold-tongue style).
TaGrid = linspace(40, 90, 24);
TbGrid = linspace(150, 330, 24);
ratioMap = nan(numel(TaGrid), numel(TbGrid));

for ia = 1:numel(TaGrid)
    for ib = 1:numel(TbGrid)
        pp = p;
        pp.T_alpha = TaGrid(ia);
        pp.T_beta = TbGrid(ib);
        pp.omega_alpha = 2*pi / pp.T_alpha;
        pp.omega_beta = 2*pi / pp.T_beta;
        pp.tspan = [0, 1800];

        [~, th] = ode45(@(tt, yy) ren2022_rhs(tt, yy, pp), pp.tspan, pp.theta0, ...
            odeset('RelTol', 1e-6, 'AbsTol', 1e-8));

        na = (th(end,1) - th(1,1)) / (2*pi);
        nb = (th(end,2) - th(1,2)) / (2*pi);
        ratioMap(ia, ib) = na / nb;
    end
end
end
