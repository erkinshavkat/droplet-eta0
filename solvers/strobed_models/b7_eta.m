function eta = b7_eta(x,y,impact,x_data,y_data,p)
    eta=0;
    for n=1:impact
        nTF=n;
        eta = eta+p.b4_prefactor/(4*pi) * (2*pi)^3 *...
                    cos(p.theta) *...
                    sqrt(pi/(p.beta1*nTF)).* ...
                    besselj(0, 2*pi .* sqrt((x - x_data(n)).^2 + (y - y_data(n)).^2 )) .* ...
                    exp(-n/p.Me -  ((x - x_data(n)).^2 + (y - y_data(n)).^2)/(4*p.beta1*nTF)) ;
    end 
end





