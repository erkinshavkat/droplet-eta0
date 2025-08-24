function [eta_b4, eta_faria] = compare_H(p)
t = p.theta/(4*pi);


H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1500, 900]); 

K_vec_scaled = p.K_vec*1.015;

H_num_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'LineWidth',2);hold on;
% H_A5_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),':',"LineWidth",2);
% H_A14_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'--',"LineWidth",2);
H_sum_ax=plot(K_vec_scaled/(2*pi),zeros(p.Nk,1),'-.',"LineWidth",2);
H_active_ax=plot(K_vec_scaled/(2*pi),zeros(p.Nk,1),'-.',"LineWidth",2);

xline(1,'k:','LineWidth',2);

xlabel('k/k_F'); ylabel('H');
legend('numerical','A13 + A14','A13*tanh + A14');

% H_num_data=zeros(p.nimpacts*p.nsteps_impact,1);
% H_A5_data=zeros(p.nimpacts*p.nsteps_impact,1);


% v = VideoWriter(sprintf('vis/A13 tanh vs no tanh theta %.2fpi.avi',p.theta/pi),'Motion JPEG AVI');
% v.FrameRate = 24;
% open(v);
L2errs=zeros(p.nimpacts*p.nsteps_impact,1);

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])

    dH_vec(:,n)=1;
    impact_ts(n)=t;
    for nn=1:p.nsteps_impact 

        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        t= t+p.dt;

        if n>p.nimpacts-3
            for impact = 1:n
                s = impact_ts(impact);
                elapsed_time = t - s;
                H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
                H_A13(:,impact) = p.H_A13(t,s,p.K_vec);
                H_A13_activation(:,impact) = p.A5_activation(elapsed_time,0.7775) *p.H_A13(t,s,p.K_vec);

                %s=3.078, m = -1.5311
            end
            sum_H_num=sum(H_vec,2);
            sum_H_formula=sum(H_A13+H_A14,2);
            sum_H_activate=sum(H_A13_activation+H_A14,2);
            



            H_num_ax.YData=sum_H_num;
            % H_A5_ax.YData=H_A5;
            % H_A14_ax.YData=H_A14;
            H_sum_ax.YData=sum_H_formula;
            H_active_ax.YData=sum_H_activate;

            L2err = norm(sum_H_num - sum_H_activate);
            L2errs((n-1)*p.nsteps_impact + nn) = L2err;

            % H_num_data((n-1)*p.nsteps_impact + nn) = H_vec;
            % H_A5_data((n-1)*p.nsteps_impact + nn) = H_A5;

            ylim([-4 4])


            title(sprintf('mem=%.2f, t=%f Tf, theta=%.2fpi', p.mem, t,p.theta/pi));

            % frame = getframe(gcf);
            % writeVideo(v, frame);
            pause(1/24); 
        end
    end
end


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





