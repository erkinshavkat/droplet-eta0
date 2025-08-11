function [eta_b4, eta_faria] = compare_H(p)
t = p.theta/(4*pi);


H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1500, 900]); 


H_num_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'LineWidth',2);hold on;
H_A5_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),':',"LineWidth",2);
H_A14_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'--',"LineWidth",2);
H_sum_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'-.',"LineWidth",2);


xline(1,'k:','LineWidth',2);

xlabel('k/k_F'); ylabel('H');
legend('numerical','A5*tanh','A14','A5*tanh+A14');

% H_num_data=zeros(p.nimpacts*p.nsteps_impact,1);
% H_A5_data=zeros(p.nimpacts*p.nsteps_impact,1);


% v = VideoWriter('vis/H long bouncer.avi','Motion JPEG AVI');
% v.FrameRate = 24;
% open(v);

L2errs=zeros(p.nimpacts*p.nsteps_impact,1);
for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
    if n==1; dH_vec(:,n)=1; end
    for nn=1:p.nsteps_impact 

        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);

        if mod(nn,1)==0
            for impact = 1
                elapsed_time = t - (impact - 1) + p.theta/(4*pi);
                H_A14(:,impact) = p.H_A14(t,p.K_vec);
                H_A5(:,impact) = p.A5_activation(elapsed_time,0.2654,1.7319) *p.H_A5(elapsed_time, p.K_vec);
            end
            sum_H_num=sum(H_vec,2);
            sum_H_formula=sum(H_A5+H_A14,2);
            H_num_ax.YData=sum_H_num;
            % H_A5_ax.YData=H_A5;
            % H_A14_ax.YData=H_A14;
            H_sum_ax.YData=sum_H_formula;

            L2err = norm(sum_H_num - sum_H_formula)/norm(sum_H_num);
            L2errs((n-1)*p.nsteps_impact + nn) = L2err;

            % H_num_data((n-1)*p.nsteps_impact + nn) = H_vec;
            % H_A5_data((n-1)*p.nsteps_impact + nn) = H_A5;

            ylim([-1 1])


            title(sprintf('mem=%.2f A5+A14, t=%f Tf, L2err=%.4f', p.mem, t, L2err));

            % frame = getframe(gcf);
            % writeVideo(v, frame);
            pause(1/24); 
        end
        t= t+p.dt;
    end


end
close all;
semilogy(L2errs)
close(v);
% close all;
% time_domain=(1:p.nimpacts*p.nsteps_impact)*p.dt;
% semilogy(time_domain,abs(H_A5_data),'LineWidth',2); hold on;
% semilogy(time_domain,abs(H_num_data),'--','LineWidth',2);
% xlabel('Time'); ylabel('H');
% legend({'H_{A13+A14}','H_{numerical}'},'Location','best');

% imagesc([0,p.Lk/(2*pi)],[0,3],H_data');
% colorbar;
% caxis([0 0.1]);
% xlabel('k/kF');
% ylabel('t/TF');
end





