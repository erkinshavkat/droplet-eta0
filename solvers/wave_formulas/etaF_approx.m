function eta = etaF_approx(x,y,t,s,p)
    %% Formerly B7, now B2/B4, approximate etaF formula from one impact
    %x: (N,1), x distance(s) from impact
    %y: (N,1), y distance(s) from impact
    %n: num of TFs since impact
    %p: params 

    eta = 2*p.b4_prefactor/(4*pi) * (p.kC)^3 *...
                cos(p.theta/2+pi/4) .*...
                cos(2*pi*t-pi/4) .*...
                sqrt(pi/(p.beta1*(t-s))).* ...
                besselj(0, p.kC .* sqrt(x.^2 + y.^2 )) .* ...
                exp(-(t-s)/p.Me -  (x.^2 + y.^2)/(4*p.beta1*(t-s))) ;
end





