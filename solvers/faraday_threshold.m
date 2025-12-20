function [is_stable,eta_max,eta_r,sqrt_energy,t_vec] = faraday_threshold(p,x)

is_stable = 1;
% Remove drop
p.drop_force{1} = @(x,y) 0;      
t=0;

% Set initial condition 
phi     = zeros(size(p.xx));
pert    = 10^(-5);
eta     = pert*( 0*exp(-5*(sqrt((p.xx-1).^2+p.yy.^2)).^2) + ...
    besselj(0,(2*pi)*sqrt((p.xx-x(1)).^2+(p.yy-x(2)).^2)) ); 

% Fourier transform (2d) of initial condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);

% storage for the energy of the wave field (energy modulo multiplicative
% constants)
sqrt_energy = zeros(1,p.nimpacts);
t_vec = zeros(1, p.nimpacts);

for n=1:p.nimpacts
   
    eta = real(ifft2(eta_hat));
    eta_r(:,:,n)=eta;
    warning('off');
    %plot_solution(eta,xi,yi,t,p)
    
    if mod(n-1,1)==0
        eta_max=max(max(abs(eta)));
        disp(['impact=',num2str(n),' eta_max=',num2str(eta_max)]);
    end
    
%     if eta_max > 10*pert
%        disp('Probably above threshold. Stopping simulation...');
%        is_stable = 0;
%        break
%     elseif eta_max < pert/10
%         disp('Probably below threshold. Stopping simulation...');
%         is_stable = 1;
%         break
% 
%     end
    
    % Define L2 norm, or the square root of the energy
    % Note: sqrt_energy should grow/decay exponentially with rate nu, where
    % nu = 0 at the Faraday threshold
    sqrt_energy(n) = sqrt(sum(abs(eta_hat(:)).^2));
    t_vec(n) = t;
        
    % Evolve wave between impacts
    [phi_hat, eta_hat] = evolve_wave_IF(phi_hat, eta_hat, t, p);
    
    t = t+p.impact_interval;

end






