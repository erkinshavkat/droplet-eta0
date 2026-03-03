function HLinf = Hinf(p)
    alpha= p.K_vec .* sqrt(p.d0* (p.K2_vec.*p.Bo + p.G));
    HLinf=  (sin(alpha)) ./ ...
        (2*alpha.* (cosh(2*p.nu0.*p.K2_vec)-cos(alpha)));

end



%sqrt(p.d0/p.Bo) *p.G /2*t