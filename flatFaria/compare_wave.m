function [x_data,y_data,t_data,eta_data,phi_hat,ui,vi] = compare_wave(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);

b1x=struct;
b1x.eta = p.eta0(:);
b1x.etaprime = zeros(size(p.eta0(:)));
b1x.xi=p.xi; b1x.yi=p.yi; b1x.ui=p.ui; b1x.vi=p.vi;
b1x.x_data = zeros(p.nimpacts,1); b1x.y_data = zeros(p.nimpacts,1); b1x.t_data = zeros(p.nimpacts,1);
b1x.eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 

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

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   
    [faria.ui, faria.vi, faria.phi_hat] = drop_impact_matt(faria.xi,faria.yi, faria.ui, faria.vi, faria.phi_hat, faria.eta_hat, p);
    [b1k.ui, b1k.vi, b1k.etaprime_hat] = b1k_impact(b1k.xi,b1k.yi, b1k.ui, b1k.vi, b1k.etaprime_hat, p);
    %[b1x.ui, b1x.vi, b1x.etaprime] = b1x_impact(b1x.xi,b1x.yi, b1x.ui, b1x.vi, b1x.eta, b1x.etaprime, p);

    %[faria.xi, faria.yi, faria.ui, faria.vi] = evolve_drops(faria.xi,faria.yi, faria.ui, faria.vi, p);
    %[b1x.xi, b1x.yi, b1x.ui, b1x.vi] = evolve_drops(b1x.xi, b1x.yi, b1x.ui, b1x.vi, p);
    

    for nn=1:p.nsteps_impact 
        [faria.phi_hat, faria.eta_hat] = evolve_wave_IF_rkstep(faria.phi_hat, faria.eta_hat, t + (nn -1)*p.dt, p); 
        [b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        %[b1x.eta, b1x.etaprime] = b1x_evolve_wave_rkstep(b1x.eta,b1x.etaprime, t + (nn -1)*p.dt, p); 


        if mod(nn,5)==0

            faria_eta=real(ifft2(faria.eta_hat));
            plot(faria_eta(:,p.Nx/2),'LineWidth',2);hold on;

            b1k.eta=real(ifft2(b1k.eta_hat));
            plot(b1k.eta(:,p.Nx/2));

            % b1x_eta_mat=reshape(b1x.eta,[p.Nx, p.Ny]);
            % plot(b1x_eta_mat(:,p.Nx/2),'LineWidth',2);

            title(['Impact ' num2str(n) ' Step ' num2str(nn)]);
            xlabel('Index'); ylabel('Value');
            ylim([-0.02,0.02])
            legend('faria','b1k')
            hold off;
            drawnow;
            pause(1/6);
        end

    end

    t = t+p.impact_interval;

end
x_data=b1x.x_data; y_data=b1x.y_data; t_data=b1x.t_data; eta_data=b1x.eta_data; ui=b1x.ui; vi=b1x.vi;
end





