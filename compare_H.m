function [eta_b4, eta_faria] = compare_H(p)

t = p.theta/(4*pi);

H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1200, 800]); 

H_full_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), 'LineWidth', 2); hold on;

% H_inf_ax = plot(p.K_vec/(2*pi), p.H_Linf(p.K_vec), 'LineStyle','--' ,'LineWidth', 2); hold on;
% H_inf_ax = plot(p.K_vec/(2*pi), p.H_Linf_approx(p.K_vec), 'LineStyle','--' ,'LineWidth', 2); hold on;

% H_taylor_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1),'LineStyle','--', 'LineWidth', 2);

% H_lin_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), ':', 'LineWidth', 2);

ylim([-1 1])
xlim([0 2])

xlabel('k/k_F'); ylabel('H');
% legend('A14 \alpha full', 'A14 \alpha quadC', 'A14 \alpha lin')
% legend('A14 accum', 'HL_{\infty}','HL_{\infty} approx')
v = VideoWriter('vis/H num.mp4','MPEG-4');
v.FrameRate = 24;
open(v);
dH_vec(:,1)=1;
xline(1, '--','Color','k');
for n=1:p.nimpacts
    disp(['Impact number: ' num2str(n)])
    % dH_vec(:,n)=1;
    impact_ts(n)=t;
    for nn=1:p.nsteps_impact
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        t= t+p.dt;

        if n>0
            % for impact = 1
            %     s = impact_ts(impact);
            %     elapsed_time = t - s;
            %     H_A14_full(:,impact)   = p.H_A14(elapsed_time,p.K_vec,p.alpha_full,p.alpha_full);
            %     H_A14_taylor(:,impact) = p.H_A14_taylor(elapsed_time,p.K_vec);
            %     % H_A14_quadC(:,impact)  = p.H_A14(elapsed_time,p.K_vec,p.alpha_quadC,p.alpha_quadC);
            %     % H_A14_lin(:,impact)    = p.H_A14(elapsed_time,p.K_vec,p.alpha_lin,p.alpha_lin);
            %     % H_A13_fixphi(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) * p.H_A13(t, s, p.K_vec, @(k) -pi/4);
            % end
            % sum_H_full   = sum(H_A14_full,2);
            % % sum_H_quadC  = sum(H_A14_quadC,2);
            % % sum_H_lin    = sum(H_A14_lin,2);
            % sum_H_taylor = sum(H_A14_taylor,2);

            H_full_ax.YData = sum(H_vec,2);
            % H_taylor_ax.YData = sum_H_taylor;
            % H_quad_ax.YData = sum_H_quadC;
            % H_lin_ax.YData  = sum_H_lin;

            frame = getframe(gcf);
            writeVideo(v, frame);

            title(sprintf('H(t,k), mem=%.2f, t=%f Tf', p.mem, t));
            % pause(1/12); 
        end
    end
end

close(v);
end





