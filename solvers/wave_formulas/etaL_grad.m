function [etaLx, etaLy] = etaL_grad(x,y,t,p)
    %x: (N,1), x distance(s) from impact
    %y: (N,1), y distance(s) from impact
    %t: time since impact
    %p: params 
    c0=2*p.d0*p.Bo + 8 *p.nu0^2;
    A=p.b4_prefactor./sqrt(2*c0*p.d0*p.Bo);
    phi=atan(2*p.nu0./sqrt(p.d0*p.Bo));
    r=sqrt(x.^2+y.^2);
    xi=- sqrt(p.d0*p.Bo).*r.^2/(t*2*c0)+sqrt(p.d0/p.Bo).*p.G/2.*t ;

    r2     = x.^2 + y.^2;
    alpha  = p.nu0 ./ (c0 .* t);
    beta   = sqrt(p.d0 .* p.Bo) ./ (2 .* c0 .* t);

    E      = exp(-alpha .* r2);
    cosXi  = cos(xi - phi);
    sinXi  = sin(xi - phi);

    % Common bracket: [-2*alpha*cos + 2*beta*sin]
    bracket = -2.*alpha .* cosXi + 2.*beta .* sinXi;

    etaLx = A./t .* E .* bracket .* x;
    etaLy = A./t .* E .* bracket .* y;
end
