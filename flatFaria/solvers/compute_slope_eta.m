function [Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p)

% Surface gradient
eta_x_hat = p.Kx.*eta_hat;
eta_y_hat = p.Ky.*eta_hat;

% Shift fourier spectrum to nearest point

ix  = find(p.xx(1,:)>xi,1); iy = find(p.yy(:,1)>yi,1);
if isempty(ix); ix=p.Nx; end
if isempty(iy); iy=p.Ny; end
shiftx = p.xx(1,ix)-xi;
shifty = p.yy(iy,1)-yi;

eta_x = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_x_hat)); % Shifts solution 
Fx = eta_x(iy,ix); 
eta_y = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_y_hat));
Fy = eta_y(iy,ix);

        
end