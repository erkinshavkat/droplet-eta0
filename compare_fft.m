function [eta_b4, eta_faria] = compare_fft(p)
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
% H_A14_full=zeros(length(p.K_vec),p.nimpacts);
% H_A14_quad=zeros(length(p.K_vec),p.nimpacts);
% H_A14_lin=zeros(length(p.K_vec),p.nimpacts); % restored
% H_A13_fixphi=zeros(length(p.K_vec),p.nimpacts);
% % H_A13_kCkf=zeros(length(p.K_vec),p.nimpacts);
% H_A13_kCform=zeros(length(p.K_vec),p.nimpacts);
% H_formula=zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1500, 900]); 

fft_domain = (0:p.Nx-1)/p.Lx;


H_full_ax = plot(fft_domain , zeros(size(fft_domain)),  'LineWidth', 2);  hold on;
% H_accumulated_ax = plot(fft_domain, zeros(size(fft_domain)), '--','LineWidth', 2);
% H_quad_ax = plot(plot_domain, zeros(size(plot_domain)), ':', 'LineWidth', 2);
% H_lin_ax = plot(plot_domain, zeros(size(plot_domain)), '--','LineWidth', 2); % restored
% H_linquad_ax = plot(plot_domain, zeros(size(plot_domain)), '-.','LineWidth', 2); % restored
xlim([0,3])
xlabel('k/k_F'); ylabel('FFT(|\eta|)');
ylim([0,0.05])
% ylim([-0.005 0.005])
% legend({'A14 \alpha(k) full','A14 \alpha(k)~sin(k^2+c)/(k^2+c)', ...
%     'A14 \alpha(k)~sin(k+k^3)/(k+k^3)','A14 H(k)~sin(k+k^3)/(k+k^3)+sin(k^2+c)/(k^2+c)'}); % updated legend
legend('A14 \alpha(k) full','A14 \alpha(k) minus 1 impact')

v = VideoWriter(sprintf('vis/b4 A14 240hz fft.avi',p.theta/pi),'Motion JPEG AVI');
v.FrameRate = 6;
open(v);

trajy=zeros(p.nimpacts,1);
dtraj=trajy(2)-trajy(1);


for n=1:p.nimpacts
    xi=0;yi=trajy(n);
    disp(['Impact number: ' num2str(n)])
    x_data(n) = xi;    y_data(n) = yi;
    [ui, vi, phi_hat] = drop_impact_matt(xi, yi, ui, vi, phi_hat, eta_hat, p);
    dH_num(:,n) = 1; 
    impact_ts(n)=t;

    for nn=1:p.nsteps_impact 
        % [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);

        if nn==p.nsteps_impact%n>p.nimpacts-3%
            for impact = 1:n
                s = impact_ts(impact);
                elapsed_time = t - s;
                alpha_quadC = @(k) p.alpha_quad(k)+sqrt(p.d0/p.Bo) *p.G /2;
                % H_A14_lin(:,impact) = p.H_A14(elapsed_time,p.K_vec,p.alpha_taylor,p.alpha_taylor); % restored
                H_A14_full(:,impact) = p.H_A14(elapsed_time,p.K_vec,p.alpha_full,p.alpha_full) ;
                % H_A14_quad(:,impact) = p.H_A14(elapsed_time,p.K_vec,alpha_quadC,alpha_quadC);
                % H_A13_fixphi(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) .* p.H_A13(t,s,p.K_vec,@(k) -pi/4);
                % H_A13_kCkf(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) .* p.H_A13(t,s,p.K_vec,@(k) p.phifunc(k,p.kf_mean));
                % H_A13_kCform(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) .* p.H_A13(t,s,p.K_vec,@(k) p.phifunc(k,p.kC_formula));
            end

            % H_A14_currentimpact = p.H_A14(t - impact_ts(end),p.K_vec,p.alpha_full,p.alpha_full); % with alpha_quad
            % eta_1impact = zeros(size(p.y));
            % eta_num_full = zeros(size(p.y));
            % eta_num_quad = zeros(size(p.y));
            % eta_num_lin = zeros(size(p.y)); % restored
            % eta_num_linquad = zeros(size(p.y)); % restored
            for i = 1:length(p.y)
                y = p.y(i);
                % eta_num_lin(i) = b4_eta(0,y, x_data,y_data,H_A14_lin,p); % restored
                % eta_1impact(i) = b4_eta(0,y, x_data,y_data,H_A14_currentimpact,p);

                eta_num_full(i) = b4_eta(0,y, x_data,y_data,H_A14_full,p);
                % eta_num_quad(i) = b4_eta(0,y, x_data,y_data,H_A14_quad,p);
                % eta_num_linquad(i) = eta_num_lin(i) + eta_num_quad(i); % restored
            end

            % FFT of eta_num_full
            eta_fft = abs(fft(eta_num_full));
            P1 = eta_fft(:,1:p.Nx/2+1);
            P1(:,2:end-1) = 2*P1(:,2:end-1);
            
            % eta_accum_fft = fftshift(abs(fft(eta_num_full-eta_1impact)));
            H_full_ax.YData = eta_fft;
            % H_accumulated_ax.YData = eta_accum_fft;

            title(sprintf('FFT(\\eta), mem=%.2f, t=%f Tf, theta=%.2fpi', p.mem, t, p.theta/pi));
            frame = getframe(gcf);
            writeVideo(v, frame);
        end
        t= t+p.dt;
    end
end

close(v);
eta_b4=[];eta_faria=[];
end





