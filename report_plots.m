clear; close all; clc

%% Load Data
fastData = load('islet_results_fast.mat');
slowData = load('islet_results_slow.mat');

fast = extract_run_metrics(fastData.t, fastData.y, fastData.p, 'Fast');
slow = extract_run_metrics(slowData.t, slowData.y, slowData.p, 'Slow');

%% Quantitative Analysis
fprintf('Fast Winding Number: %.3e\n', fast.slope)
fprintf('Slow Winding Number: %.3f\n', slow.slope)

%% Figure 1: Winding number plots using terminal window
fig1 = figure('Color', 'w', 'Units','inches','Position',[1 1 4.5 5]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Fast regime
ax1 = nexttile;
hold(ax1,'on');

plot(ax1, fast.theta_b_mean_unwrap, fast.theta_a_mean_unwrap, ...
    'k', 'LineWidth', 0.85);

plot(ax1, fast.theta_b_mean_unwrap(fast.fit_idx), ...
    fast.theta_a_mean_unwrap(fast.fit_idx), ...
    'k', 'LineWidth', 1.5);

xfit = linspace(min(fast.theta_b_mean_unwrap(fast.fit_idx)), ...
                max(fast.theta_b_mean_unwrap(fast.fit_idx)), 200);
yfit = polyval(fast.pf, xfit);
plot(ax1, xfit, yfit, 'r--', 'LineWidth', 2.5);

xlabel(ax1, '\theta_\beta');
ylabel(ax1, '\theta_\alpha');
title(ax1, 'Fast Regime');
text(ax1, 0.05, 0.95, ...
    sprintf('slope = %.4f', fast.slope), ...
    'Units','normalized', 'VerticalAlignment','top', ...
    'FontSize',10, 'BackgroundColor','w');
box(ax1,'off');

% Slow regime
ax2 = nexttile;
hold(ax2,'on');

plot(ax2, slow.theta_b_mean_unwrap, slow.theta_a_mean_unwrap, ...
    'k', 'LineWidth', 0.85);

plot(ax2, slow.theta_b_mean_unwrap(slow.fit_idx), ...
    slow.theta_a_mean_unwrap(slow.fit_idx), ...
    'k', 'LineWidth', 1.5);

xfit = linspace(min(slow.theta_b_mean_unwrap(slow.fit_idx)), ...
                max(slow.theta_b_mean_unwrap(slow.fit_idx)), 200);
yfit = polyval(slow.pf, xfit);
plot(ax2, xfit, yfit, 'r--', 'LineWidth', 2.5);

xlabel(ax2, '\theta_\beta');
ylabel(ax2, '\theta_\alpha');
title(ax2, 'Slow Regime');
text(ax2, 0.05, 0.95, ...
    sprintf('slope = %.4f', slow.slope), ...
    'Units','normalized', 'VerticalAlignment','top', ...
    'FontSize',10, 'BackgroundColor','w');
box(ax2,'off');

sgtitle('Winding Numbers, Terminal 25% of Simulation');

exportgraphics(fig1,'Figures/fig_winding_terminal.png','Resolution',300);

%% Figure 2: Mean Phase Oscillations

fig2 = figure('Color','w','Units','inches','Position',[1 1 4.5 5.0]);
t2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

xlims = [min(fast.t) max(fast.t)];
ylims = [0 2*pi];

% Fast regime
ax1 = nexttile;
hold(ax1,'on');
plot(ax1, fast.t, fast.theta_b_mean_wrap, 'r', 'LineWidth', 1.2);
plot(ax1, fast.t, fast.theta_a_mean_wrap, 'g', 'LineWidth', 1.2);
xlim(ax1, xlims);
ylim(ax1, ylims);
ylabel(ax1, 'Average \theta');
title(ax1, 'Fast Regime');
legend(ax1, '\theta_\beta', '\theta_\alpha', 'Location', 'northwest');
box(ax1,'off');
set(ax1, 'FontSize', 9, ...
    'YTick', [0 pi 2*pi], ...
    'YTickLabel', {'0', '\pi', '2\pi'});

% Slow regime
ax2 = nexttile;
hold(ax2,'on');
plot(ax2, slow.t, slow.theta_b_mean_wrap, 'r', 'LineWidth', 1.2);
plot(ax2, slow.t, slow.theta_a_mean_wrap, 'g', 'LineWidth', 1.2);
xlim(ax2, xlims);
ylim(ax2, ylims);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Average \theta');
title(ax2, 'Slow Regime');
box(ax2,'off');
set(ax2, 'FontSize', 9, ...
    'YTick', [0 pi 2*pi], ...
    'YTickLabel', {'0', '\pi', '2\pi'});

exportgraphics(fig2,'Figures/fig_phase_oscillations_wrapped.png','Resolution',300);


%% Figure 3a: Glucagon concentration snapshots, fast regime
run = fast;
Npts = 4;
t_pts = linspace(run.t(1), run.t(end), Npts);

fig3a = figure('Color','w','Units','inches','Position',[1 1 4.5 4]);
t3a = tiledlayout(2,2,'TileSpacing','compact');

clims = [0, max(fast.Gmax, slow.Gmax)];

for i = 1:Npts
    ax = nexttile;
    hold(ax,'on');

    [~, j] = min(abs(run.t - t_pts(i)));
    G_i = reshape(run.G_all(j,:), run.p.Ny, run.p.Nx);

    pcolor(ax, run.p.X, run.p.Y, G_i);
    shading(ax, 'interp');
    view(ax, 2);
    axis(ax, 'equal');
    axis(ax, 'tight');
    clim(ax, clims);

    scatter(ax, run.p.xa, run.p.ya, 2, 'g', 'filled');
    scatter(ax, run.p.xb, run.p.yb, 2, 'r', 'filled');

    title(ax, sprintf('t = %.0f s', run.t(j)));

    % Only left column gets y-label
    if ismember(i, [1 3])
        ylabel(ax, 'y (\mum)');
    else
        ylabel(ax, '');
        ax.YTickLabel = [];
    end

    % Only bottom row gets x-label
    if ismember(i, [3 4])
        xlabel(ax, 'x (\mum)');
    else
        xlabel(ax, '');
        ax.XTickLabel = [];
    end

    box(ax,'off');
    set(ax,'FontSize',8);

    % Panel label on first tile
    if i == 1
        text(ax, -0.3, 1.45, '(A)', ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold', ...
            'FontSize',11, ...
            'Clipping','off');
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Glucagon concentration';

hA = scatter(nan, nan, 18, 'g', 'filled');
hB = scatter(nan, nan, 18, 'r', 'filled');

lgd = legend([hA, hB], {'\alpha cells', '\beta cells'}, ...
    'Orientation','horizontal', ...
    'Box','off');
lgd.Layout.Tile = 'north';

sgtitle('Fast Regime');

exportgraphics(fig3a,'Figures/fig_3a_glucagon_fast.png','Resolution',300);

%% Figure 3b: Glucagon concentration snapshots, slow regime
run = slow;
Npts = 4;
t_pts = linspace(run.t(1), run.t(end), Npts);

fig3b = figure('Color','w','Units','inches','Position',[1 1 4.5 4]);
t3b = tiledlayout(2,2,'TileSpacing','compact');

clims = [0, max(fast.Gmax, slow.Gmax)];

for i = 1:Npts
    ax = nexttile;
    hold(ax,'on');

    [~, j] = min(abs(run.t - t_pts(i)));
    G_i = reshape(run.G_all(j,:), run.p.Ny, run.p.Nx);

    pcolor(ax, run.p.X, run.p.Y, G_i);
    shading(ax, 'interp');
    view(ax, 2);
    axis(ax, 'equal');
    axis(ax, 'tight');
    clim(ax, clims);

    scatter(ax, run.p.xa, run.p.ya, 2, 'g', 'filled');
    scatter(ax, run.p.xb, run.p.yb, 2, 'r', 'filled');

    title(ax, sprintf('t = %.0f s', run.t(j)));

    % Only left column gets y-label
    if ismember(i, [1 3])
        ylabel(ax, 'y (\mum)');
    else
        ylabel(ax, '');
        ax.YTickLabel = [];
    end

    % Only bottom row gets x-label
    if ismember(i, [3 4])
        xlabel(ax, 'x (\mum)');
    else
        xlabel(ax, '');
        ax.XTickLabel = [];
    end

    box(ax,'off');
    set(ax,'FontSize',8);

    % Panel label on first tile
    if i == 1
        text(ax, -0.3, 1.45, '(B)', ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold', ...
            'FontSize',11, ...
            'Clipping','off');
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Glucagon concentration';

hA = scatter(nan, nan, 18, 'g', 'filled');
hB = scatter(nan, nan, 18, 'r', 'filled');

lgd = legend([hA, hB], {'\alpha cells', '\beta cells'}, ...
    'Orientation','horizontal', ...
    'Box','off');
lgd.Layout.Tile = 'north';

sgtitle('Slow Regime');

exportgraphics(fig3b,'Figures/fig_3b_glucagon_slow.png','Resolution',300);

%% Helper function: Data Extraction from Simulation output
function run = extract_run_metrics(t, y, p, label)

    run = struct();
    run.label = label;
    run.t = t;
    run.p = p;

    % Phase blocks from state vector
    theta_a_raw = y(:, 1:p.Na);
    theta_b_raw = y(:, p.Na + (1:p.Nb));

    % Wrapped phases for display
    theta_a_wrap = mod(theta_a_raw, 2*pi);
    theta_b_wrap = mod(theta_b_raw, 2*pi);

    % Unwrapped phases for trend / winding analysis
    theta_a_unwrap = unwrap(theta_a_raw, [], 1);
    theta_b_unwrap = unwrap(theta_b_raw, [], 1);

    % Mean phase traces
    theta_a_mean_wrap = mean(theta_a_wrap, 2);
    theta_b_mean_wrap = mean(theta_b_wrap, 2);

    theta_a_mean_unwrap = mean(theta_a_unwrap, 2);
    theta_b_mean_unwrap = mean(theta_b_unwrap, 2);

    % Wrapped phase difference
    dtheta = theta_b_mean_wrap - theta_a_mean_wrap;

    % Wrap phase difference to [-pi, pi]
    dtheta = atan2(sin(dtheta), cos(dtheta));

    % Winding-number fit from unwrapped means at terminal window index
    % To ignore initial transient behavior, focus on the last 25% of the
    % simulation so we calculate from settled behavior.
    frac_terminal = 0.25;
    n = numel(t);
    idx0 = max(1, floor((1 - frac_terminal)*n));
    fit_idx = idx0:n;
    pf = polyfit(theta_b_mean_unwrap(fit_idx), theta_a_mean_unwrap(fit_idx), 1);
    slope = pf(1);
    intercept = pf(2);

    % Glucagon block
    G_all = y(:, p.Na + p.Nb + 1:end);

    % Store outputs
    run.theta_a_raw = theta_a_raw;
    run.theta_b_raw = theta_b_raw;

    run.theta_a_wrap = theta_a_wrap;
    run.theta_b_wrap = theta_b_wrap;

    run.theta_a_unwrap = theta_a_unwrap;
    run.theta_b_unwrap = theta_b_unwrap;

    run.theta_a_mean_wrap = theta_a_mean_wrap;
    run.theta_b_mean_wrap = theta_b_mean_wrap;

    run.theta_a_mean_unwrap = theta_a_mean_unwrap;
    run.theta_b_mean_unwrap = theta_b_mean_unwrap;

    run.dtheta = dtheta;

    run.fit_idx = fit_idx;
    run.fit_t0 = t(idx0);
    run.frac_terminal = frac_terminal;

    run.slope = slope;
    run.intercept = intercept;
    run.pf = pf;

    run.G_all = G_all;
    run.Gmax = max(G_all(:));
end