function etaF = etaF_bouncer(x,y,p)
    %% Formerly B10, now C1, approximate etaF formula from a bouncer
    %x: (N,1), x distance(s) from impact
    %y: (N,1), y distance(s) from impact
    %p: params 
    etaF=p.b4_prefactor/(4) * (2*pi)^3 *...
                    cos(p.theta) *...
                    sqrt(p.Me/p.beta1).* ...
                    besselj(0, 2*pi .* sqrt(x.^2 + y.^2 )) .* ...
                    exp(-  sqrt((x.^2 + y.^2)/(p.Me*p.beta1))) ;
end





