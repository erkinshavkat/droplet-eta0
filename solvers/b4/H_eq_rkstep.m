function [Hfunc,Hprime] = H_eq_rkstep(Hfunc,Hprime, t_in, p)
% one step of the Runge-Kutta method for the damped mathieu



t = t_in;    


[rhs_hfunc_1, rhs_hprime_1] = H_rhs(Hfunc, ...
                                      Hprime, ...
                                      t, p);

[rhs_hfunc_2, rhs_hprime_2] = H_rhs(Hfunc      + p.dt/2 * rhs_hfunc_1, ...
                                    Hprime + p.dt/2 * rhs_hprime_1, ...
                                    t+p.dt/2, p);

[rhs_hfunc_3, rhs_hprime_3] = H_rhs(Hfunc      + p.dt/2 * rhs_hfunc_2, ...
                                    Hprime + p.dt/2 * rhs_hprime_2, ...
                                    t+p.dt/2, p);


[rhs_hfunc_4, rhs_hprime_4] = H_rhs(Hfunc      + p.dt * rhs_hfunc_3, ...
                                    Hprime + p.dt * rhs_hprime_3, ...
                                    t+p.dt, p);


Hfunc      = Hfunc      + p.dt / 6 * (rhs_hfunc_1 + 2*rhs_hfunc_2 + 2*rhs_hfunc_3 + rhs_hfunc_4);
Hprime = Hprime + p.dt / 6 * (rhs_hprime_1 + 2*rhs_hprime_2 + 2*rhs_hprime_3 + rhs_hprime_4);


end

     
