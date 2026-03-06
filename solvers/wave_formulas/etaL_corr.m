function etaL = etaL_corr(x,y,t,p)
    %% B11, approximate etaL formula from one impact
    %x: (N,1), x distance(s) from impact
    %y: (N,1), y distance(s) from impact
    %t: time since impact
    %p: params 
    %this ends up being the exact same integral as etaL for a gaussian spatial impact
    %only difference is R^2/2 is replaced by a/b
    b = p.G/2*sqrt(p.d0/p.Bo);
    a = sqrt(p.d0*p.Bo);
    R=sqrt(2*a/b);
    A=p.b4_prefactor/a;
    qt = 4*t^2*(a^2+4.*p.nu0^2)+R^4+8*p.nu0*R^2*t;
    r2=x.^2+y.^2;
    xi = t.*(b-a*r2./(qt));
    etaL = -A * exp(-r2*(R^2+4*p.nu0*t)/(2*qt)).* ...
        ((R^2+4*p.nu0*t).*sin(xi) + 2*t.*a.*cos(xi))./qt;
end



%sqrt(p.d0/p.Bo) *p.G /2*t