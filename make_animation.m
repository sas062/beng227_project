% MP4 & GIF of glucagon concentration over time

clc; clear; close all;
load('islet_results.mat');

figure;
g_all = y(:, p.Na + p.Nb + 1:end);
cmax = max(g_all(:));

v = VideoWriter('./Figures/glucagon_animation.mp4', 'MPEG-4');
v.FrameRate = 15;
open(v);

for k = 1:length(t)
    Gk = reshape(y(k, p.Na + p.Nb + 1:end), p.Ny, p.Nx);

    imagesc(p.x, p.y, Gk);
    axis xy equal tight;
    clim([0 cmax]);
    colorbar;
    hold on;
    plot(p.xa, p.ya, 'wo', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    plot(p.xb, p.yb, 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    hold off;

    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(sprintf('Glucagon Concentration, t = %.2f', t(k)));

    drawnow;
    writeVideo(v, getframe(gcf));
end

close(v);

filename = './Figures/glucagon_animation.gif';

figure;
g_all = y(:, p.Na + p.Nb + 1:end);
cmax = max(g_all(:));

for k = 1:length(t)
    Gk = reshape(y(k, p.Na + p.Nb + 1:end), p.Ny, p.Nx);

    imagesc(p.x, p.y, Gk);
    axis xy equal tight;
    clim([0 cmax]);
    colorbar;
    hold on;
    plot(p.xa, p.ya, 'wo', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    plot(p.xb, p.yb, 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    hold off;

    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(sprintf('Glucagon Concentration, t = %.2f', t(k)));

    drawnow;

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.08);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.08);
    end
end
