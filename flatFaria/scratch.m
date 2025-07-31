clear all; close all; clc;
load("results/b4vfaria.mat")
eta_b4_1 = eta_b4;
eta_faria_1 = eta_faria;

load("results/b4vfaria_2x.mat")
eta_b4_2 = eta_b4;
eta_faria_2 = eta_faria;

nframes = min([size(eta_b4_1,2), size(eta_faria_1,2), size(eta_b4_2,2), size(eta_faria_2,2)]);

v = VideoWriter('compare_all_eta.avi','Motion JPEG AVI');
v.FrameRate = 3;
open(v);
domain=linspace(-8,8,128);
for frame = 1:nframes
    clf;
    plot(domain,eta_b4_1(:,frame),'LineWidth',2); hold on;
    plot(domain,eta_b4_2(:,frame),'LineWidth',2);
    plot(domain,eta_faria_1(:,frame),'--','LineWidth',2);
    ylim([-0.01 0.01]);
    xlim([-1,1])
    legend({'b4 ','b4\_refine', ...
            'faria'}, ...
            'Location','best');
    title(['Frame ' num2str(frame)]);
    frame_img = getframe(gcf);
    writeVideo(v, frame_img);
end

close(v);