clc; clear; close all;
load('islet_results_fast.mat');

%% Plot beta-cell and alpha-cell phases

% Extract all phase histories
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));

theta_a_all = mod(theta_a_all, 2*pi);
theta_b_all = mod(theta_b_all, 2*pi);

figure(1);
plot(t,theta_a_all);
xlabel("time (s)")
ylabel("\theta_\alpha")

figure(2);
plot(t,theta_b_all);
xlabel("time (s)")
ylabel("\theta_\beta")

%% Plot mean beta-cell and mean alpha-cell traces
% Extract all phase histories
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));
theta_a_all = mod(theta_a_all, 2*pi);
theta_b_all = mod(theta_b_all, 2*pi);

theta_b_mean = mean(theta_b_all,2);
theta_a_mean = mean(theta_a_all,2);

figure(3);
plot(t,theta_a_mean,'g');
xlabel("time (s)")
ylabel("Average \theta_\alpha")

figure(4);
plot(t,theta_b_mean,'r');
xlabel("time (s)")
ylabel("Average \theta_\beta")

figure(5);
hold on;
plot(t,theta_b_mean,'r');
plot(t,theta_a_mean,'g');
ylabel("Average \theta")
legend("\theta_\beta","\theta_\alpha")

%% Plot glucagon concentrations at six timepoints
Npts = 6;
t_pts = linspace(0,max(t),Npts);
fig3 = figure('Color','w', 'Position', [100 100 900 500]);
hold on;
for i = 1:Npts
    subplot(2,3,i)
    hold on;
    j = find(t == t_pts(i));
    G_i = reshape(y(j, p.Na + p.Nb + 1:end), p.Ny, p.Nx);
    G_all = y(:,p.Na + p.Nb + 1:end);

    pcolor(p.X,p.Y,G_i);
    shading interp;
    view(2);
    axis equal tight;
    cb2 = colorbar;
    clim([0, max(G_all(:))])
    title(sprintf('Glucagon Concentration\n t = %d', t_pts(i)));
    scatter(p.xa,p.ya,5*ones(p.Na,1),'filled','g')
    scatter(p.xb,p.yb,5*ones(p.Nb,1),'filled','r')
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    %legend([hA2, hB2], {'Alpha cells', 'Beta cells'}, 'Location', 'best');
end

%% Plot Global Winding number
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));
% theta_a_all = mod(theta_a_all, 2*pi);
% theta_b_all = mod(theta_b_all, 2*pi);

theta_b_mean = mean(theta_b_all,2);
theta_a_mean = mean(theta_a_all,2);

figure(7)
hold on;
plot(theta_b_mean,theta_a_mean,'k');

% Best-fit line
pf = polyfit(theta_b_mean, theta_a_mean, 1);
xfit = linspace(min(theta_b_mean), max(theta_b_mean), 200);
yfit = polyval(pf, xfit);
slope = pf(1);
intercept = pf(2);

fprintf('Best-fit slope (theta_a vs theta_b): %.6f\n', slope);
fprintf('Best-fit intercept: %.6f\n', intercept);
plot(xfit, yfit, 'r--', 'LineWidth', 2);
ylabel('\theta_\alpha')
xlabel('\theta_\beta')


% Plot phase difference

theta_a_all = mod(theta_a_all, 2*pi);
theta_b_all = mod(theta_b_all, 2*pi);

theta_b_mean = mean(theta_b_all,2);
theta_a_mean = mean(theta_a_all,2);
dtheta = theta_b_mean-theta_a_mean;
figure(8);
plot(t,dtheta,'k');
ylabel('\Delta\theta')
xlabel("time (s)")


%% Plot Global Synchronicity between alpha and beta cells
theta_a_all = y(:, 1:p.Na);
theta_b_all = y(:, p.Na + (1:p.Nb));

theta_b_mean = mod(mean(theta_b_all,2),2*pi);
theta_a_mean = mod(mean(theta_a_all,2),2*pi);



R = synchronicity([theta_a_mean theta_b_mean]);

figure(9);
plot(t,R);
xlabel("time (s)");
ylabel("Kuramoto Order Parameter (R)")

%% B-B coupling 
load('islet_results.mat');
if p.diffusion == 0
    
    theta_b_all = y(:, p.Na + (1:p.Nb));
    
    theta_b_mean = mod(mean(theta_b_all,2),2*pi);
    
    
    
    
    R = synchronicity([theta_b_all]);
    
    figure();
    plot(t,R);
    xlabel("time (s)");
    ylabel("Kuramoto Order Parameter (R)")
end
function R = synchronicity(theta)
    [~,N] = size(theta);
    R_tmp = 0;
    for i = 1:N
        R_tmp = R_tmp + exp(1i*theta(:,i));
    end
    R = abs(R_tmp)/N;
end
