% Perplexity code...
% MP4 & lightweight GIF of glucagon concentration over time

clc; clear; close all;
load('islet_results.mat');

outDir = './Figures';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

g_all = y(:, p.Na + p.Nb + 1:end);
cmax = max(g_all(:));

nFramesTotal = length(t);

%% -------------------------
%  MP4: more detailed version
%  Target duration ~ 20 s
%  -------------------------
mp4_target_duration = 20;                 % seconds
mp4_fps = 20;                             % playback fps
mp4_target_frames = round(mp4_target_duration * mp4_fps);
mp4_stride = max(1, floor(nFramesTotal / mp4_target_frames));
mp4_idx = 1:mp4_stride:nFramesTotal;

fig1 = figure('Color','w', 'Position', [100 100 900 700]);
v = VideoWriter(fullfile(outDir, 'glucagon_animation.mp4'), 'MPEG-4');
v.FrameRate = mp4_fps;
open(v);

for k = mp4_idx
    Gk = reshape(y(k, p.Na + p.Nb + 1:end), p.Ny, p.Nx);

    clf(fig1);
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
    title(sprintf('Glucagon Concentration, t = %.1f s', t(k)));

    drawnow;
    writeVideo(v, getframe(fig1));
end

close(v);

%% -------------------------
%  GIF: lightweight preview
%  Target duration ~ 12–15 s
%  Fewer frames, smaller image, fewer colors
%  -------------------------
gif_target_duration = 12;                 % seconds
gif_delay = 0.12;                         % seconds per GIF frame
gif_target_frames = round(gif_target_duration / gif_delay);
gif_stride = max(1, floor(nFramesTotal / gif_target_frames));
gif_idx = 1:gif_stride:nFramesTotal;

gif_filename = fullfile(outDir, 'glucagon_animation_preview.gif');

fig2 = figure('Color','w', 'Position', [150 150 560 420]);

for ii = 1:length(gif_idx)
    k = gif_idx(ii);
    Gk = reshape(y(k, p.Na + p.Nb + 1:end), p.Ny, p.Nx);

    clf(fig2);
    imagesc(p.x, p.y, Gk);
    axis xy equal tight;
    clim([0 cmax]);
    colorbar;
    hold on;
    plot(p.xa, p.ya, 'wo', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
    plot(p.xb, p.yb, 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    hold off;

    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title(sprintf('Glucagon, t = %.0f s', t(k)));

    drawnow;

    frame = getframe(fig2);
    im = frame2im(frame);

    % Shrink GIF frame to reduce file size
    im = imresize(im, 0.55);

    % Reduce color count for smaller GIF
    [imind, cm] = rgb2ind(im, 64);

    if ii == 1
        imwrite(imind, cm, gif_filename, 'gif', ...
            'Loopcount', inf, 'DelayTime', gif_delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', ...
            'WriteMode', 'append', 'DelayTime', gif_delay);
    end
end
