function [ui, vi, etaprime] = b1x_impact(xi,yi, ui, vi, eta,etaprime, p)

eta_mat=reshape(eta,[p.Nx, p.Ny]);
eta_hat=fft2(eta_mat);
% DESCRIBES THE INSTANTANEOUS INTERACTION BETWEEN THE DROP AND THE SURFACE
[Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p);

% Drop speed immediately after impact
ui = - p.G.*Fx./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*ui;
vi = - p.G.*Fy./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*vi;
% Wave immediately after impact


etaprime=etaprime +  p.Pvec(xi,yi);

end