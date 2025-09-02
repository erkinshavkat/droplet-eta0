function [eta_b4, eta_faria] = compare_wave(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);




phi = p.phi0; 
eta = p.eta0; 
phi_hat = fft2(phi); 
eta_hat = fft2(eta);
xi = p.xi; yi = p.yi; ui = p.ui; vi = p.vi;

x_data= zeros(p.nimpacts,1);
y_data= zeros(p.nimpacts,1);


H_num=zeros(length(p.K_vec),p.nimpacts);
dH_num= zeros(length(p.K_vec),p.nimpacts);
H_A14=zeros(length(p.K_vec),p.nimpacts);
H_A13=zeros(length(p.K_vec),p.nimpacts);
H_A13_activate=zeros(length(p.K_vec),p.nimpacts);


H_formula=zeros(length(p.K_vec),p.nimpacts);


%faria_ax=plot(p.y,zeros(p.Nx,1),'LineWidth',2);hold on;

plotting=false;

fig = figure('Position', [0, 0, 1500, 900]); 


plot_range = 8;

plot_domain=-plot_range:p.hy:plot_range;
H_num_ax=plot(plot_domain,0*plot_domain,'LineWidth',2);hold on;
% H_A5_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),':',"LineWidth",2);
% H_A14_ax=plot(p.K_vec/(2*pi),zeros(p.Nk,1),'--',"LineWidth",2);
% H_test_ax=plot(plot_domain,0*plot_domain,'-.',"LineWidth",2);
H_active_ax=plot(plot_domain,0*plot_domain,'-.',"LineWidth",2);

xlim([-plot_range,plot_range])

legend('numerical','A14','A13*tanh + A14');
v = VideoWriter('vis/b4 A14only bouncer.avi','Motion JPEG AVI');
v.FrameRate = 24; % 6 frames per second, adjust as needed
open(v);

num_waveheight=[];
formula_waveheight=[];
activate_waveheight=[];
t_vec=[];

trajy=zeros(p.nimpacts,1);%linspace(-1,1,p.nimpacts);
dtraj=trajy(2)-trajy(1);

K_vec_scaled = p.K_vec;
K3_vec_scaled = K_vec_scaled.^3;
for n=1:p.nimpacts
    xi=0;yi=trajy(n);
    disp(['Impact number: ' num2str(n)])
    x_data(n) = xi;    y_data(n) = yi;

    [ui, vi, phi_hat] = drop_impact_matt(xi, yi, ui, vi, phi_hat, eta_hat, p);

    dH_num(:,n) = 1; 
    impact_ts(n)=t;



    for nn=1:p.nsteps_impact 

        %[phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        %[b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);

        
        if n > p.nimpacts-5
            for impact = 1:n
                s = impact_ts(impact);
                elapsed_time = t - s;
                H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
                H_A13(:,impact) = p.H_A13(t,s,p.K_vec);
            end

            b4_eta_compute = @(x, y, impact, H) p.b4_prefactor * sum(K3_vec_scaled.* H(:,impact) .* ...
            besselj(0, K_vec_scaled .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
            
            
            eta_centerline = @(y,H) sum(arrayfun(@(impact) b4_eta_compute(p.x(p.Nx/2), y, impact, H), 1:n));

            % H_splice = H_num;
            % splice_split=0.8;
            % mask = p.K_vec <= splice_split*2*pi;
            % H_splice(mask,:)=H_A5_activate(mask,:)+H_A14(mask,:);

            
            eta_num = arrayfun(@(y)eta_centerline(y,H_num), plot_domain);
            eta_sum = arrayfun(@(y)eta_centerline(y,H_A14+H_A13), plot_domain);
            % eta_sum= arrayfun(@(y)eta_centerline(y,H_A13+H_A14), plot_domain);

            %faria_ax.YData=faria_plot;
            H_num_ax.YData=eta_num;
            H_active_ax.YData=eta_sum;
            % H_test_ax.YData=eta_sum;

            ylim([-0.02 0.02])


            %title(sprintf('mem=%.2f, t=%f Tf, splice [formula|%.2f kF|numerical]', p.mem, t,splice_split));
            title(sprintf('mem=%.2f, t=%f Tf, theta=%.2fpi', p.mem, t,p.theta/pi));
            frame = getframe(gcf);
            writeVideo(v, frame);
            % pause(1/24); 
        end
        t= t+p.dt;


        % if n>27
        %     next_x = xi ; next_y=yi+dtraj;
        %     for impact = 1:n
        %         s = impact_ts(impact);
        %         elapsed_time = t - s;
        %         H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
        %         H_A13(:,impact) = p.H_A13(t,s,p.K_vec);
        %         H_A13_activate(:,impact) = p.A5_activation(elapsed_time,0.7775).*p.H_A13(t,s,p.K_vec);
        %     end
        %     b4_eta_compute = @(x, y, impact, H) p.b4_prefactor * sum(p.K3_vec .* H(:,impact) .* besselj(0, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        %     eta_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_eta_compute(x, y, impact, H), 1:n));
        %     num_waveheight=([num_waveheight,eta_centerline(next_x,next_y,H_num)]);
        %     formula_waveheight=([formula_waveheight,eta_centerline(next_x,next_y,H_A13+H_A14)]);
        %     activate_waveheight=([activate_waveheight,eta_centerline(next_x,next_y,H_A13_activate+H_A14)]);
        %     t_vec=[t_vec,t];
        % end

    end

end

% plot(t_vec, num_waveheight,"LineWidth",2);hold on
% title(sprintf('eta difference at (0,0), theta=%.2fpi', p.theta/pi));
% plot(t_vec, formula_waveheight,"LineWidth",2);
% plot(t_vec, activate_waveheight,"LineWidth",2);
% for i=28:p.nimpacts
%     xline( p.theta/(4*pi) + i ,"LineStyle",':','LineWidth',2)
% end
% legend('num','no tanh','tanh')
% ylim([0 0.012])

% saveas(gcf,sprintf('vis/Wave error 00 theta %.2fpi.png',p.theta/pi))
close(v);
% close all;

% imagesc([0,p.Lk/(2*pi)],[0,3],H_data');
% colorbar;
% caxis([0 0.1]);
% xlabel('k/kF');
% ylabel('t/TF');
end





