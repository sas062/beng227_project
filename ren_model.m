p.w_a = 2*pi / 40;
p.w_b = 2*pi / 360;
p.K_ba = -2 * pi / 10;
p.K_ab = 2*pi / 10;
p.tau = 0.5;

y0 = [1; 0; 0]; % th_a0, th_b0, G0

tspan = [0, 1000];

[t,y] = ode23(@(t,y) rhs(t,y,p), tspan, y0);

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

figure;
hold on
plot(t, mod(y(:,1),2*pi), 'g', 'LineWidth', 1.5)
plot(t, mod(y(:,2),2*pi), 'r', 'LineWidth', 1.5)
plot(t, y(:,3), 'w', 'LineWidth', 1.5)
%ylim([0 2*pi]);
hold off

xlabel('Time')
legend('\theta_a', '\theta_b', 'G');


