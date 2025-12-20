function HL = HL_formula(t,k,p)
%%HL given by A11
%inputs:
% t: time elapsed since impact (equivalent to t-s in HF, here it's simplified bc we only have t-s dependency)
% k: k domain, can be (Nk,1)
% p: params struct
alpha= k.* sqrt(p.d0 *(k.^2*p.Bo + p.G));
% alpha_quad =sqrt(d0*Bo).*k.^2;
% alpha_lin = sqrt(d0*G).*k;
% alpha_linquad = sqrt(d0*G).*k + sqrt(d0*Bo).*k.^2;
% alpha_taylor = sqrt(d0*G).*k + sqrt(d0/G)*Bo/2.*k.^3;
% alpha_quadC = alpha_quad(k)+sqrt(d0/Bo) *G /2;
% alpha_linquadC = alpha_lin(k) + alpha_quadC(k);

HL = exp(-2*p.nu0*k.^2.*t).*sin(alpha.*t)./alpha;
end

