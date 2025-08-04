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


H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

faria_ax=plot(p.y,zeros(p.Nx,1),'LineWidth',2);hold on;
b4_ax=plot(p.y,zeros(p.Nx,1),'--',"LineWidth",2);
xlabel('Index'); ylabel('Value');
xlim([-3,3])
legend('faria','b4')

% v = VideoWriter('moving b4.avi','Motion JPEG AVI');
% v.FrameRate = 12; % 6 frames per second, adjust as needed
% open(v);

eta_b4=zeros(p.Nx,p.nimpacts*p.nsteps_impact); eta_faria = zeros(p.Nx,p.nimpacts*p.nsteps_impact);
H_data = zeros(length(p.K_vec),p.nimpacts*p.nsteps_impact);
for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
    yi=(n-1)*0.5-2;
    x_data(n) = xi;    y_data(n) = yi;

    [ui, vi, phi_hat] = drop_impact_matt(xi, yi, ui, vi, phi_hat, eta_hat, p);
    dH_vec(:,n) = 1; 

    for nn=1:p.nsteps_impact 

        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        %[b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);


        if mod(nn,1)==0

            b4_eta_compute = @(x, y, impact) p.b4_prefactor * sum(p.K3_vec .* H_vec(:,impact) .* besselj(0, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
            b4_eta_center = @(y) sum(arrayfun(@(impact) b4_eta_compute(p.x(p.Nx/2), y, impact), 1:n));
            b4_eta = arrayfun(b4_eta_center, p.y);

            faria_eta = real(ifft2(eta_hat));
            faria_plot = faria_eta(:,p.Nx/2);

            b4_eta_plot=b4_eta;
            faria_ax.YData=faria_plot;
            b4_ax.YData=b4_eta_plot;
            
            eta_b4(:,(n-1)*p.nsteps_impact + nn ) = b4_eta_plot;
            eta_faria(:,(n-1)*p.nsteps_impact + nn ) = faria_plot;
            ylim([-0.005 0.005])


            title(['Impact ' num2str(n) ' Step ' num2str(nn)]);

            % frame = getframe(gcf);
            % writeVideo(v, frame);
            pause(1/48); 
        end
        t= t+p.dt;
    end


end


close(v);
close all;

% imagesc([0,p.Lk/(2*pi)],[0,3],H_data');
% colorbar;
% caxis([0 0.1]);
% xlabel('k/kF');
% ylabel('t/TF');
end





