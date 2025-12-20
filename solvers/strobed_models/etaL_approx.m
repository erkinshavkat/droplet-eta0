function etaL = etaL_approx(x,y,xp,yp,p,t)
    c0=2*p.d0*p.Bo + 8 *p.nu0^2;
    A=p.b4_prefactor./sqrt(2*c0*p.d0*p.Bo);
    phi=atan(2*p.nu0./sqrt(p.d0*p.Bo));
    r=sqrt((x-xp).^2+(y-yp).^2);
    xi=- sqrt(p.d0*p.Bo).*r.^2/(t*2*c0)+sqrt(p.d0/p.Bo).*p.G/2.*t ;


    etaL = A * exp(-p.nu0.*r.^2/(c0*t)) .* cos(xi - phi)/t;
end



%sqrt(p.d0/p.Bo) *p.G /2*t