function [HLinf,HLinf_approx] = Hinftau(tau,p)
    alpha_full= p.K_vec .* sqrt(p.d0 *(p.K2_vec.*p.Bo + p.G));
    alpha_quad = sqrt(p.d0*p.Bo).*p.K_vec.^2;
    alpha_quadC = alpha_quad+sqrt(p.d0/p.Bo) *p.G /2;
    alpha_lin =  sqrt(p.d0*p.G).*p.K_vec;

    alpha=alpha_full;
    alpha_recip=alpha_full;
    k=p.K_vec;
    k2=p.K2_vec;
    HLinf = exp(-2*p.nu0*k2.*tau) .* ...
        (sin(alpha.*(1-tau)) + exp(2*p.nu0.*k2).*sin(alpha.*tau)) ./ ...
        (2*alpha_recip.* (cosh(2*p.nu0.*k2)-cos(alpha)));

    alpha=alpha_full;
    alpha_recip=alpha_quad;
    k=p.K_vec;
    k2=p.K2_vec;


    HLinf_approx = sin(alpha).*exp(-2*p.nu0*k2) ./ alpha .*(1/2 + cos(alpha)./(p.d0*p.G*k2));
 
end



%sqrt(p.d0/p.Bo) *p.G /2*t