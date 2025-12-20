function [x_data,y_data,t_data,eta_data,ui,vi] = b1x_trajectory(p)
%% Set initial condition 
% <<< MATT >>> I've changed the initial impact time to account for the
% change in gravitational acceleration. Impacts happen at times t_n = n +
% t_0 (in dimensionless variables), with t_0 = theta/(4*pi) (theta = impact
% phase). The new impact times are still stored in the vector t_data.
% t = 0; % <<< OLD VALUE >>>
t = p.theta/(4*pi);

eta = p.eta0(:); 
etaprime = zeros(size(eta));

xi = p.xi; yi=p.yi; ui=p.ui; vi=p.vi; 

% Create vector for position
x_data = zeros(p.nimpacts,1); y_data = x_data; t_data = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 

%[ui, vi, etaprime] = impact_b1_xspace(xi,yi, ui, vi, etaprime, p);

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   
    eta_max=max(max(abs(eta)));

    if eta_max > 1
        disp('Probably above threshold. Stopping simulation...');
        break
    end
    disp(xi);

    % Store position (and possibly wavefield)
    x_data(n,:) = xi; y_data(n,:) = yi; t_data(n) = t; eta_data(:,:,n) = gather(reshape(eta,[p.Nx,p.Ny]));
    
    % Drop impact
    
    [ui, vi, etaprime] = b1x_impact(xi,yi, ui, vi,eta, etaprime, p);
    
    % Evolve drops between impacts
    [xi, yi, ui, vi] = evolve_drops(xi, yi, ui, vi, p);
    
    % Evolve wave between impacts
    for nn=1:p.nsteps_impact 
        [eta, etaprime] = b1x_evolve_wave_rkstep(eta,etaprime, t + (nn -1)*p.dt, p); 

    end
    %disp(['wave',num2str(max(max(ksi_hat)))])

    t = t+p.impact_interval;

end
end





