function [etaLx, etaLy] = etaL_corr_grad(x,y,t,p)
    %% vibecoded etaL correction gradient, from etaL_corr formula
    % constants are probably all over the place
    b  = p.G/2 * sqrt(p.d0/p.Bo);
    a  = sqrt(p.d0 * p.Bo);
    R  = sqrt(2*a/b);
    A  = p.b4_prefactor / a;

    qt = 4*t^2*(a^2 + 4.*p.nu0^2) + R^4 + 8*p.nu0*R^2*t;
    r2 = x.^2 + y.^2;

    xi   = t .* (b - a.*r2./qt);
    gam  = (R^2 + 4*p.nu0*t);          % scalar shorthand
    F    = exp(-r2 .* gam ./ (2*qt));

    sinXi = sin(xi);
    cosXi = cos(xi);

    S    = gam .* sinXi + 2*t*a .* cosXi;

    % Derivatives with respect to r2
    dxi_dr2 = -t*a / qt;                              % scalar
    dF_dr2  = F .* (-gam / (2*qt));                   % array
    dS_dr2  = dxi_dr2 .* (gam .* cosXi - 2*t*a .* sinXi);  % array

    % d/dr2 [ F * S ] = dF_dr2 * S + F * dS_dr2
    % d/dx  [ F * S ] = d/dr2[F*S] * 2x

    dFSdr2 = dF_dr2 .* S + F .* dS_dr2;

    etaLx = -A/qt .* dFSdr2 .* (2.*x);
    etaLy = -A/qt .* dFSdr2 .* (2.*y);
end



%sqrt(p.d0/p.Bo) *p.G /2*t