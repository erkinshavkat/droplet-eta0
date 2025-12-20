function HLinf = Hinftau(tau,k,p)
    %% C4 formula for Hinf
    k2=p.K2_vec;
    alpha_full= k .* sqrt(p.d0 *(k2.*p.Bo + p.G));
    % alpha_quad = sqrt(p.d0*p.Bo).*k2;
    % alpha_quadC = alpha_quad+sqrt(p.d0/p.Bo) *p.G /2;
    % alpha_lin =  sqrt(p.d0*p.G).*k;

    alpha=alpha_full;
    alpha_recip=alpha_full;

    HLinf = exp(-2*p.nu0*k2.*tau) .* ...
        (sin(alpha.*(1-tau)) + exp(2*p.nu0.*k2).*sin(alpha.*tau)) ./ ...
        (2*alpha_recip.* (cosh(2*p.nu0.*k2)-cos(alpha)));
 
end

