function [ui, vi, etaprime_hat] = b1k_impact(xi,yi, ui, vi, etaprime_hat, p)

%{
% DESCRIBES THE INSTANTANEOUS INTERACTION BETWEEN THE DROP AND THE SURFACE
[Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p);

% Drop speed immediately after impact
ui = - p.G.*Fx./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*ui;
vi = - p.G.*Fy./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*vi;
% Wave immediately after impact
%}

etaprime_hat=etaprime_hat +  p.d0*p.M(1)/(p.hx*p.hy)*p.G*p.K2_deriv.*exp(-p.Kx*(p.Lx/2+xi)-p.Ky*(p.Ly/2+yi));

end