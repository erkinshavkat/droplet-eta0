function [eta_b4, eta_faria] = compare_wave(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);
disp(p.nsteps_impact)
phi = p.phi0; 
eta = p.eta0; 
phi_hat = fft2(phi); 
eta_hat = fft2(eta);
xi = p.xi; yi = p.yi; ui = p.ui; vi = p.vi;

x_data= zeros(p.nimpacts,1);
y_data= zeros(p.nimpacts,1);

H_num=zeros(length(p.K_vec),p.nimpacts);
dH_num= zeros(length(p.K_vec),p.nimpacts);
H_A14_full=zeros(length(p.K_vec),p.nimpacts);
% H_A14_quad=zeros(length(p.K_vec),p.nimpacts);
% H_A14_approx=zeros(length(p.K_vec),p.nimpacts);
% H_A14_taylor=zeros(length(p.K_vec),p.nimpacts);
% H_A13_fixphi=zeros(length(p.K_vec),p.nimpacts);
H_A13_kCkf=zeros(length(p.K_vec),p.nimpacts);
% H_A13_kCform=zeros(length(p.K_vec),p.nimpacts);
% H_formula=zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1200, 900]); 

plot_range = 2;
plot_domain = -plot_range:p.hy:plot_range;


% H_num_ax = plot(plot_domain, zeros(size(plot_domain)),  'LineWidth', 2);  hold on;
H_A14_full_ax = plot(plot_domain, zeros(size(plot_domain)), 'LineWidth', 2); hold on;
alphak_approx_ax = plot(plot_domain, zeros(size(plot_domain)), '--', 'LineWidth', 2); hold on;
% H_A13_kCkf_ax = plot(plot_domain, zeros(size(plot_domain)), '--', 'LineWidth', 2); hold on;
% sum_ax= plot(plot_domain, zeros(size(plot_domain)), ':', 'LineWidth', 2); hold on;
% H_A14_taylor_ax = plot(plot_domain, zeros(size(plot_domain)), '-.', 'LineWidth', 2); hold on;
% H_A14_inf_ax = plot(plot_domain, eta_A14_inf, 'LineWidth', 2); hold on; % static plot
% H_A14_asymp_ax = plot(plot_domain, zeros(size(plot_domain)), 'LineWidth', 2); hold on; % static plot
legend({'\alpha\sim k','Approx','Sum'});
xlim([-plot_range,plot_range])
ylim([-0.005 0.005])

v = VideoWriter(sprintf('vis/test.avi',p.theta/pi),'Motion JPEG AVI');
v.FrameRate = 6;
open(v);

trajy=zeros(p.nimpacts,1);

for n=1:p.nimpacts
    xi=0;yi=trajy(n);
    disp(['Impact number: ' num2str(n)])
    x_data(n) = xi;    y_data(n) = yi;
    % dH_num(:,n) = 1; 
    impact_ts(n)=t;
    for nn=1:p.nsteps_impact
        % [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);

        if true
            eta_L_formula = zeros(size(plot_domain));
            for impact = 1
                s = impact_ts(impact);
                elapsed_time = t - s;
                H_A14_full(:,impact)   = p.H_A14(elapsed_time,p.K_vec,p.alpha_full,p.alpha_full);
                % H_A13_kCkf(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) .* p.H_A13(t,s,p.K_vec,@(k) p.phifunc(k,p.kf_mean));
                % H_A14_taylor(:,impact) = p.H_A14(elapsed_time,p.K_vec,p.alpha_taylor,p.alpha_taylor);
                eta_L_formula = eta_L_formula + etaL_approx(0,plot_domain,0,0,p,elapsed_time);

            end

            % eta_A14_full=zeros(p.Nx,1);
            % eta_A13=zeros(p.Nx,1);
            for i=1:length(plot_domain)
                y = plot_domain(i);

                eta_A14_full(i)   = eta_num(0,y, x_data,y_data,H_A14_full,p);
                % eta_A13(i)       = eta_num(0,y, x_data,y_data,H_A13_kCkf,p);
                % % eta_A14_lin(i)    = eta_num(0,y, x_data,y_data,H_A14_lin,p);
                % eta_A14_taylor(i) = eta_num(0,y, x_data,y_data,H_A14_taylor,p);
            end
            H_A14_full_ax.YData   = eta_A14_full;
            alphak_approx_ax.YData = eta_L_formula;
            % H_A13_kCkf_ax.YData  = eta_A13;
            % sum_ax.YData = eta_A14_full + eta_A13;
            % H_A14_quad_ax.YData   = eta_A14_quad;
            % % H_A14_lin_ax.YData    = eta_A14_lin;
            % % H_A14_taylor_ax.YData = eta_A14_taylor;
            % H_A14_asymp_ax.YData  = eta_L_formula;
            pause(1/24);
            % frame = getframe(gcf);
            % title(sprintf('\\eta_L, t=%.2f Tf', t));
            % writeVideo(v, frame);
        end
        t= t+p.dt;

    end
end
close(v);
% eta_b4=[];eta_faria=[];
end





