function [d_eta, d_etaprime] = b1x_eta_rhs(eta, etaprime, t, p)

    
    d_eta=etaprime;
    d_etaprime = 4*p.nu0* p.Lap2d*etaprime ...
     -(4*p.nu0^2 + p.d0*p.Bo)*p.Lap2d*p.Lap2d*eta  ...
          +p.d0*p.g(t)*p.Lap2d*eta;
end