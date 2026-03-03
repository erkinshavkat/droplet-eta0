function HL = HL_formula(varargin)

t=varargin{1};k=varargin{2};p=varargin{3};
alpha= k.* sqrt(p.d0 *(k.^2*p.Bo + p.G));

if nargin>3
    alphatype=varargin{4};
    switch alphatype
        case 'quad'
            alpha = sqrt(p.d0*p.Bo).*k.^2;
        case 'quadC'
            alpha = sqrt(p.d0*p.Bo).*k.^2+sqrt(p.d0/p.Bo)*p.G/2;
        case 'lin'
            alpha = sqrt(p.d0*p.G).*k;
    end

end
alpha_recip=alpha;
if nargin==5
    alphareciptype=varargin{5};
    switch alphareciptype
        case 'quad'
            alpha_recip = sqrt(p.d0*p.Bo).*k.^2;
        case 'quadC'
            alpha_recip = sqrt(p.d0*p.Bo).*k.^2+sqrt(p.d0/p.Bo)*p.G/2;
        case 'lin'
            alpha_recip = sqrt(p.d0*p.G).*k;
    end
end
% alpha_quad = @(k) sqrt(d0*Bo).*k.^2;
% alpha_lin = @(k) sqrt(d0*G).*k;
% alpha_linquad = @(k) sqrt(d0*G).*k + sqrt(d0*Bo).*k.^2;
% alpha_taylor = @(k) sqrt(d0*G).*k + sqrt(d0/G)*Bo/2.*k.^3;
% alpha_quadC = @(k) alpha_quad(k)+sqrt(d0/Bo) *G /2;
% alpha_linquadC = @(k) alpha_lin(k) + alpha_quadC(k);

HL = exp(-2*p.nu0*k.^2.*t).*sin(alpha.*t)./alpha_recip;
end

