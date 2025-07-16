function [eta,etaprime] = b1x_evolve_wave_rkstep(eta,etaprime, t_in, p)
% Evolves wave field in F-space through nsteps of size dt

t = t_in;    
    % stage 1
[d_eta_1, d_etaprime_1] = b1x_eta_rhs(eta, ...
                                      etaprime, ...
                                      t, p);

[d_eta_2, d_etaprime_2] = b1x_eta_rhs(eta      + p.dt/2 * d_eta_1, ...
                                      etaprime + p.dt/2 * d_etaprime_1, ...
                                      t+p.dt/2, p);

[d_eta_3, d_etaprime_3] = b1x_eta_rhs(eta      + p.dt/2 * d_eta_2, ...
                                      etaprime + p.dt/2 * d_etaprime_2, ...
                                      t+p.dt/2, p);


[d_eta_4, d_etaprime_4] = b1x_eta_rhs(eta      + p.dt * d_eta_3, ...
                                      etaprime + p.dt * d_etaprime_3, ...
                                      t+p.dt, p);



eta      = eta      + p.dt/6* (d_eta_1      + 2* d_eta_2      + 2* d_eta_3      + d_eta_4 );
etaprime = etaprime + p.dt/6* (d_etaprime_1 + 2* d_etaprime_2 + 2* d_etaprime_3 + d_etaprime_4 );



    
end

     
