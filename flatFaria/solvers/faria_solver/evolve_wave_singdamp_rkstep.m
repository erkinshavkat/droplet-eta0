function [phi_hat, eta_hat] = evolve_wave_singdamp_rkstep(phi_hat, eta_hat,t_in,p)
% Evolves wave field in F-space through nsteps of size dt

t = t_in;

% we use a modified Runge-Kutta method for applying the dissipation (given
% by the matrix D) analytically. Ask MD for notes.

    % stage 1
[rhs1_1, rhs2_1] = compute_rhs_full_IF(phi_hat, eta_hat, t, p);


% stage 2
[rhs1_2, rhs2_2] = compute_rhs_full_IF(phi_hat + p.dt/2.*rhs1_1, ...
    eta_hat + p.dt/2.*rhs2_1, t+p.dt/2, p); 

% stage 3
[rhs1_3, rhs2_3] = compute_rhs_full_IF(phi_hat + p.dt/2.*rhs1_2, ...
    eta_hat + p.dt/2.*rhs2_2, t+p.dt/2, p); 

% stage 4
[rhs1_4, rhs2_4] = compute_rhs_full_IF( (phi_hat + p.dt.*rhs1_3), ...
    ( eta_hat + p.dt.*rhs2_3), t+p.dt, p);

% RK step
phi_hat = p.D.^(2) .* (phi_hat + p.dt/6 * (rhs1_1 + 2*rhs1_2 + 2*rhs1_3 + rhs1_4));
eta_hat = p.D.^(2) .* (eta_hat + p.dt/6 * (rhs2_1 + 2*rhs2_2 + 2*rhs2_3 + rhs2_4));
   
end

     
