function [x_data,y_data,eta_data,eta_intermediate] = trajectory_IF_matt(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);

phi = p.phi0;
eta = p.eta0; 
xi = p.xi; yi=p.yi; ui=p.ui; vi=p.vi; 

% Create vector for position
x_data = zeros(p.nimpacts,1); y_data = x_data; t_data = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 
eta_intermediate= zeros(p.Nx,p.Ny,p.nsteps_impact); 

% Fourier transform (2d) of initial condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   
    eta = real(ifft2(eta_hat));   
    eta_max=max(max(abs(eta)));
    if eta_max > 1
        disp('Probably above threshold. Stopping simulation...');
        break
    end   
    % Store position (and possibly wavefield)
    x_data(n,:) = xi; y_data(n,:) = yi; t_data(n) = t;
    eta_data(:,:,n ) = eta;
    % Drop impact
    % <<< MATT >>> use the updated function
    [ui, vi, phi_hat] = drop_impact_matt(xi,yi, ui, vi, phi_hat, eta_hat, p);
    
    % Evolve drops between impacts
    %[xi, yi, ui, vi] = evolve_drops(xi, yi, ui, vi, p);
    % Evolve wave between impacts
    for nn=1:p.nsteps_impact 
        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t + (nn -1)*p.dt, p); 
        if n>p.nimpacts-1
            eta_intermediate(:,:,nn ) = gather(real(ifft2(eta_hat)));
        end
    end
    t = t+p.impact_interval;

end
end





