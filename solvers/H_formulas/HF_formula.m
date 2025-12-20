function HF = HF_formula(varargin)
%%HF given by A8, with choices of phi, default is constant -pi/4
% Inputs:
% t: time of evaluation
% s: time of impact
% k: k domain, can be (Nk,1)
% p: params struct
% phi: optional, constant if nothing, 'linear' for linear approximation, 'cotan' for full cotan formula
phi=-pi/4;

t=varargin{1};s=varargin{2};k=varargin{3};p=varargin{4};
if nargin>4
    phitype=varargin{5};
    switch phitype
        case 'linear'
            phi = -pi/4+1/2*(sqrt(p.beta1/(p.nu0*(p.kC)^2)).*(k-p.kC)+ 2*p.nu0*p.kC^2*(p.mem-1) );
        case 'cotan'
            phi = mod(-1/2*acot(sqrt(p.beta1/(p.nu0*(p.kC)^2)).*(k-p.kC)+ 2*p.nu0*p.kC^2*(p.mem- 1) ),-pi/2);
    end
end


beta_func= - 1/(p.Me)-p.beta1.*(k-p.kC).^2;



HF= 2/(4*pi) *exp(beta_func.*(t-s)).*cos(2*pi*t + phi).*cos(2*pi*s-phi);

end