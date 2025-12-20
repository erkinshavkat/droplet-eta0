function etaL = etaL_approx(x,y,t,p)
    %% B11, approximate etaL formula from one impact
    %x: (N,1), x distance(s) from impact
    %y: (N,1), y distance(s) from impact
    %t: time since impact
    %p: params 
    c0=2*p.d0*p.Bo + 8 *p.nu0^2;
    A=p.b4_prefactor./sqrt(2*c0*p.d0*p.Bo);
    phi=atan(2*p.nu0./sqrt(p.d0*p.Bo));
    r=sqrt(x.^2+y.^2);
    xi=- sqrt(p.d0*p.Bo).*r.^2/(t*2*c0)+sqrt(p.d0/p.Bo).*p.G/2.*t ;


    etaL = A * exp(-p.nu0.*r.^2/(c0*t)) .* cos(xi - phi)/t;
end



%sqrt(p.d0/p.Bo) *p.G /2*t