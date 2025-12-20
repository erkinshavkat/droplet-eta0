function [etaL1, etaL2] = etaLinf(x,y,xp,yp,tau,p)

    c0=2*p.d0*p.Bo + 8 *(p.nu0)^2;
    A=p.b4_prefactor./sqrt(4*p.d0*p.Bo);
    phi=atan(2*p.nu0./sqrt(p.d0*p.Bo));
    r=sqrt((x-xp).^2+(y-yp).^2);
    xi=- sqrt(p.d0*p.Bo).*r.^2/(tau*2*c0)+sqrt(p.d0/p.Bo).*p.G/2.*tau ;

    mu = 8*p.nu0^2+(1+tau)^2 + p.d0*p.Bo*(1-tau)^2;

    phitilde = atan(2*p.nu0./sqrt(p.d0*p.Bo)*(1+tau)/(1-tau));
    xitilde =- sqrt(p.d0*p.Bo).*(1-tau).*r.^2/(2*mu)+sqrt(p.d0/p.Bo).*p.G/2.*(1-tau) ;

    etaL1 = A .* (sqrt(1/(2*c0)).*exp(-p.nu0.*r.^2./(c0*tau)) .* cos(xi - phi)./tau-...
        exp(-p.nu0.*(1-tau).*r.^2./mu) .* cos(xitilde - phitilde)./sqrt(1-tau));

    etaL2 = A .* (sqrt(1/(2*c0)).*exp(-p.nu0.*r.^2./(c0*tau)) .* cos(xi - phi)./tau);
end



%sqrt(p.d0/p.Bo) *p.G /2*t