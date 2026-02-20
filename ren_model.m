clc; clear all; clf; close all;

p.w_a = 2*pi / 40;
p.w_b = 2*pi / 360;
p.tau = 0.5;

% K_ab, K_ba
% K = [
%     0.17, -0.40
%     0.80, -0.75
%     0.88, -0.41
%     0.13, -0.04
% ];

K = [
    0.17,  -0.66   % red    
    0.52,  -0.66   % green
    0.33,  -0.66   % blue
    0.084, -0.37   % purple
];

K = K .* (0.2*pi);

titles = ["Slow","Fast","Mixed","Weakly Coupled"];

% th_a0, th_b0, G0
Y0 = [
    pi, 0, 0
    pi/2,   0, 0
    3*pi/4,      0, 0
    3*pi/2,     0, 0
]; 

fig = figure;
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
ax = gobjects(4,1);

tspan = [0, 360];

for i = 1:size(K,1)
    p.K_ab = K(i,1);
    p.K_ba = K(i,2);
    y0 = Y0(i,:).';
  
    [t,y] = ode23(@(t,y) rhs(t,y,p), tspan, y0);

    % figure(i)
    % hold on
    % plot(t, mod(y(:,1),2*pi), 'g', 'LineWidth', 1.5)
    % plot(t, mod(y(:,2),2*pi), 'r', 'LineWidth', 1.5)
    % %plot(t, y(:,3), 'w', 'LineWidth', 1.5)
    % ylim([0 2*pi]);    
    % xlabel('Time')
    % ylabel('\theta');
    % title(titles(i));
    % legend('\theta_a', '\theta_b');
    % hold off

    ax(i) = nexttile(tl); hold(ax(i),'on')
    plot(ax(i), t, mod(y(:,1),2*pi), 'g', 'LineWidth', 1.5)
    plot(ax(i), t, mod(y(:,2),2*pi), 'r', 'LineWidth', 1.5)

    ylim(ax(i), [0 2*pi]);
    ylabel(ax(i), '\theta');
    title(ax(i), titles(i));
    if i == 1
        legend(ax(i), '\theta_a', '\theta_b', 'Location','best');
    end
    if i < 4
        ax(i).XTickLabel = []; 
    else
        xlabel(ax(i), 'Time');
    end
    hold(ax(i),'off')

end

function dthdt = rhs(~,y,p)

th_a = y(1);
th_b = y(2);
G = y(3);

fs = @(th) (1-cos(th)) / 2;
fra = @(th) sin(th);
frb = @(th) 1;

dth_a = p.w_a + p.K_ba * fs(th_b) * fra(th_a);
dth_b = p.w_b + p.K_ab * G * frb(th_b);
dG = fs(th_a) - p.tau * G;

dthdt = [dth_a; dth_b; dG];
end



