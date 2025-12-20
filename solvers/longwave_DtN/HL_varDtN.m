function HL = HL_varDtN(t,k,p,DtN)

alpha= sqrt(DtN.*(k.^2.*p.Bo + p.G));
% alpha_quad = @(k) sqrt(d0*Bo).*k.^2;
% alpha_lin = @(k) sqrt(d0*G).*k;
% alpha_linquad = @(k) sqrt(d0*G).*k + sqrt(d0*Bo).*k.^2;
% alpha_taylor = @(k) sqrt(d0*G).*k + sqrt(d0/G)*Bo/2.*k.^3;
% alpha_quadC = @(k) alpha_quad(k)+sqrt(d0/Bo) *G /2;
% alpha_linquadC = @(k) alpha_lin(k) + alpha_quadC(k);


HL = exp(-2*p.nu0*k.^2.*t).*sin(alpha.*t)./alpha;
end

