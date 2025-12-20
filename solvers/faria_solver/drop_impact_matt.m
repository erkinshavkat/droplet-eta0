function [ui, vi, phi_hat] = drop_impact_matt(xi,yi, ui, vi, phi_hat, eta_hat, p)
% DESCRIBES THE INSTANTANEOUS INTERACTION BETWEEN THE DROP AND THE SURFACE
[Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p);

% Drop speed immediately after impact
ui = - p.G.*Fx./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*ui;
vi = - p.G.*Fy./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*vi;
% Wave immediately after impact

% <<< MATT >>> I've changed the jump condition to account for the minus
% sign in front of the pressure term
% phi_hat = phi_hat + p.M(1)*p.G/(p.hx*p.hy).*exp(-p.Kx*(p.Lx/2+xi(1))-p.Ky*(p.Ly/2+yi(1))); % <<< OLD VERSION >>>
phi_hat = phi_hat - p.M(1)*p.G/(p.hx*p.hy).*exp(-p.Kx*(p.Lx/2+xi(1))-p.Ky*(p.Ly/2+yi(1)));

end