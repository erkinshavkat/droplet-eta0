function [x_data,y_data,t_data,eta_data,phi_hat,ui,vi] = compare_wave(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);

faria=struct;
faria.phi = p.phi0; faria.eta = p.eta0; 
faria.phi_hat = fft2(faria.phi); faria.eta_hat = fft2(faria.eta);
faria.xi=p.xi; faria.yi=p.yi; faria.ui=p.ui; faria.vi=p.vi;
faria.x_data = zeros(p.nimpacts,1); faria.y_data = zeros(p.nimpacts,1); faria.t_data = zeros(p.nimpacts,1);
faria.eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 



b1k=struct;
b1k.xi=p.xi; b1k.yi=p.yi; b1k.ui=p.ui; b1k.vi=p.vi;
b1k.etaprime = zeros(size(p.xx));
b1k.eta = p.eta0; 
b1k.eta_hat=fft2(b1k.eta);
b1k.etaprime_hat=fft(b1k.etaprime);
b1k.x_data = zeros(p.nimpacts,1); b1k.y_data = zeros(p.nimpacts,1); b1k.t_data = zeros(p.nimpacts,1);
b1k.eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 


[faria.ui, faria.vi, faria.phi_hat] = drop_impact_matt(faria.xi,faria.yi, faria.ui, faria.vi, faria.phi_hat, faria.eta_hat, p);
[b1k.ui, b1k.vi, b1k.etaprime_hat] = b1k_impact(b1k.xi,b1k.yi, b1k.ui, b1k.vi, b1k.etaprime_hat, p);



prefactor=  -p.d0*p.M(1)*p.G/(2*pi);

H_vec=zeros(p.Nx*p.Ny,1);
dH_vec= ones(p.Nx*p.Ny,1);


b4_waveheight=zeros(p.nimpacts*p.nsteps_impact);
faria_waveheight=zeros(p.nimpacts*p.nsteps_impact);
disp(p.x(p.Nx/2+1))
for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   


    for nn=1:p.nsteps_impact 
        [faria.phi_hat, faria.eta_hat] = evolve_wave_IF_rkstep(faria.phi_hat, faria.eta_hat, t + (nn -1)*p.dt, p); 
        %[b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t + (nn -1)*p.dt, p);


        % faria_eta=real(ifft2(faria.eta_hat));
        % b4_eta_00=prefactor *sum(p.K3_vec.* H_vec.*besselj(0,p.K_vec.* sqrt((0-faria.xi).^2 +(0-faria.yi).^2 )));
        
        % b4_waveheight((n-1)*p.nimpacts + nn)=b4_eta_00;
        % faria_waveheight((n-1)*p.nimpacts + nn))
        if mod(nn,5)==0

            b4_eta_plot  = zeros(size(p.y)); 
            for index=1:length(p.y)
                b4_eta_plot(index) = prefactor *sum(p.K3_vec.* H_vec.*besselj(0,p.K_vec.* sqrt((p.x(p.Nx/2)-faria.xi).^2 +(p.y(index)-faria.yi).^2 )));
            end
            faria_eta=real(ifft2(faria.eta_hat));
            faria_plot=faria_eta(:,p.Nx/2);

            % plot(p.y,faria_plot,'LineWidth',2);hold on;
            % plot(p.y,b4_eta_plot);
            % ylim([-0.02,0.02])
            plot(p.y,faria_plot/max(abs(faria_plot)),'LineWidth',2);hold on;
            plot(p.y,b4_eta_plot/max(abs(b4_eta_plot)));
            ylim([-1,1])
            title(['Impact ' num2str(n) ' Step ' num2str(nn)]);
            xlabel('Index'); ylabel('Value');

            xlim([-3,3])
            legend('faria','b4')
            hold off;
            drawnow;
            pause(1/6);
        end

    end

    t = t+p.impact_interval;

end
x_data=b1x.x_data; y_data=b1x.y_data; t_data=b1x.t_data; eta_data=b1x.eta_data; ui=b1x.ui; vi=b1x.vi;
end





