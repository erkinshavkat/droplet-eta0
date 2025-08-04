function [eta_b4, eta_faria] = compare_H(p)
t = p.theta/(4*pi);


H_vec=zeros(length(p.K_vec),1);
dH_vec= ones(length(p.K_vec),1);

H_num_ax=plot(p.K_vec/(2*pi),H_vec,'LineWidth',2);hold on;
H_A5_ax=plot(p.K_vec/(2*pi),H_vec,'--',"LineWidth",2);
xlabel('Index'); ylabel('Value');

legend('faria','b4')

% v = VideoWriter('moving b4.avi','Motion JPEG AVI');
% v.FrameRate = 12; % 6 frames per second, adjust as needed
% open(v);
for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])

    for nn=1:p.nsteps_impact 

        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);


        if mod(nn,1)==0
            H_A5= arrayfun (@(k) p.H_A5(t *p.TF,k) +p.Hlong2(t*p.TF, k), p.K_vec /(2*pi)*p.kf_mean);

            H_num_ax.YData=H_vec/max(H_vec);
            H_A5_ax.YData=H_A5/max(H_A5);


            ylim([-0.5 0.5])


            title(['Impact ' num2str(n) ' Step ' num2str(nn)]);

            % frame = getframe(gcf);
            % writeVideo(v, frame);
            pause(1/24); 
        end
        t= t+p.dt;
    end


end


% close(v);
close all;

% imagesc([0,p.Lk/(2*pi)],[0,3],H_data');
% colorbar;
% caxis([0 0.1]);
% xlabel('k/kF');
% ylabel('t/TF');
end





