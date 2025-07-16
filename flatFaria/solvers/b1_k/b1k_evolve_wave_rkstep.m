function [eta_hat,etaprime_hat] = b1k_evolve_wave_rkstep(eta_hat,etaprime_hat, t_in, p)
% Evolves wave field in F-space through nsteps of size dt



t = t_in;    
Deta_hat = eta_hat; Detaprime_hat = p.kappa.*eta_hat + etaprime_hat;


[rhs_eta_1, rhs_etaprime_1] = b1k_eta_rhs(Deta_hat, ...
                                      Detaprime_hat, ...
                                      t, p);

[rhs_eta_2, rhs_etaprime_2] = b1k_eta_rhs(Deta_hat      + p.dt/2 * rhs_eta_1, ...
                                          Detaprime_hat + p.dt/2 * rhs_etaprime_1, ...
                                          t+p.dt/2, p);

[rhs_eta_3, rhs_etaprime_3] = b1k_eta_rhs(Deta_hat      + p.dt/2 * rhs_eta_2, ...
                                          Detaprime_hat + p.dt/2 * rhs_etaprime_2, ...
                                          t+p.dt/2, p);


[rhs_eta_4, rhs_etaprime_4] = b1k_eta_rhs(Deta_hat      + p.dt * rhs_eta_3, ...
                                          Detaprime_hat + p.dt * rhs_etaprime_3, ...
                                          t+p.dt, p);


RK_eta_term=p.dt / 6 * (rhs_eta_1 + 2*rhs_eta_2 + 2*rhs_eta_3 + rhs_eta_4);
RK_etaprime_term=p.dt / 6 * (rhs_etaprime_1 + 2*rhs_etaprime_2 + 2*rhs_etaprime_3 + rhs_etaprime_4);

Dinv_eta_hat = RK_eta_term; 
Dinv_etaprime_hat = -p.kappa.*RK_eta_term + RK_etaprime_term;

eta_hat      = p.D.^2 .* (eta_hat      + Dinv_eta_hat);
etaprime_hat = p.D.^2 .* (etaprime_hat + Dinv_etaprime_hat);


end

     
