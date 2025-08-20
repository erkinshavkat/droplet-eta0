function p = setup_IF_matt(Gam,H,eta0,Nx,Lx,Nk,kmin,kmax)
% Sets most of the parameters for the problem
% Input: 
%   Nx          -------- Number of points in x
%   Ny          -------- Number of points in x
%   Lx          -------- Size of domain in x
%   Ly          -------- Size of domain in y
%   Gam         -------- Amplitude of the shaking
%   dt_desired  -------- Desired time step
%% Set parameters
mem =0.98;
Gam = mem*Gam;
Ny = Nx; 
Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/16;


%% Grid, variable coefficient, and initial data:
Nxy = Nx*Ny;
hx = Lx/Nx; hy = Ly/Ny; Lx = Lx; Ly = Ly;
x = hx*(0:Nx-1)-Lx/2; y = hy*(0:Ny-1)-Ly/2;
[xx,yy] = meshgrid(x,y);
[xxx,yyy] = meshgrid(-Lx/2:0.02:Lx/2-hx,-Ly/2:0.02:Ly/2-hy);

%% Problem parameters in SI units (typically do not change for a given experiment)
nu        = 2*10^(-5);        % Kinematic viscosity (m2/s)
                              % Effective viscosity
nu        = 0.865*nu;       % 70Hz - Shallow h=1.5 mm (og. 0.865)
rho       = 949;              % Density (kg/m3);
sig       = 0.0206;           % Surface tension (N/m);
omega0     = 70*2*pi;         % Angular frequency in 1/s
g0         = 9.81;             % Gravity in m/s

%% Dispersion relation for inviscid surface waves
% <<< MATT >>> changed to include the tanh within the square root
% dispEuler = @(k,h)sqrt( g0 .* k + sig./rho*k.^3 ).* tanh(k*h); % <<< OLD VALUE >>>
dispEuler = @(k,h)sqrt( g0 .* k + sig./rho*k.^3 ).* sqrt(tanh(k*h)); 

%% Finds wave number of subharmonic mode in deep and shallow regions 
Gamma_neutral = @(k,h) sqrt(4./(g0*k*tanh(k*h)).^2.*( (dispEuler(k,h).^2+(2*nu*k.^2).^2 ...
                                    -(omega0/2)^2).^2 + omega0.^2*(2*nu*k.^2).^2));

[kf_mean   ,Gamma_max_mean] = fminsearch(@(k)Gamma_neutral(k,H),1300);
kf_mean = kf_mean*1.02;
b0=(tanh(kf_mean*H) / kf_mean);

% for i=1:10
% kf_mean = sqrt(rho*omega0^2/(2*sqrt(b0)*(rho*g0*sqrt(b0)+sqrt(sig*rho*omega0^2+b0*rho^2*g0^2))));
% b0=(tanh(kf_mean*H) / kf_mean);
% end


kf_deep=kf_mean*ones(Nx,Ny);
Gamma_max_deep=Gamma_max_mean*ones(Nx,Ny);
h_top_grid=H*ones(Nx,Ny);
options = optimset('Display','off');

%% Computes "effective depth" of model
lambdaf = 2*pi/kf_mean;
d_deep  = tanh(kf_deep.*h_top_grid) ./ kf_deep;

%% Drop parameters (dimensional)
drop_radius = (0.745/2)*10^(-3); 

% <<< MATT >>> I've changed the value of the impact phase to account for
% the change in definition of vibrational acceleration below. The new value
% of theta might need to be tweaked, but it should be about the value that
% I've listed below. I expect that we will have 1.35 <= theta/pi <= 1.55.
% The closer theta gets to pi, the faster the droplet walks. The further
% theta gets from pi, the slower the droplet walks. I used theta = 1.43 *
% pi in the corral paper, for example.

%theta = 0.266*2*pi; %(0.042 xF/TF) at Me 0.99 % <<< OLD VALUE >>>
%theta = (1 + 2 * 0.248) * pi; % <<< try this value first >>>
theta = 1.5*pi; % <<< matt's testing value >>>

drop_density  = 949;              % Density of drop (kg/m3)
drop_mass = 4/3*pi*drop_radius^3*drop_density; % mass of drop (kg);

%% Air viscosity
mu_air      = 1.8*10^(-5);               % Viscosity of air [kg / (m s)]

%% Choice of scales
TF          = 4*pi/omega0;      % Chosen time scale (Faraday period)
xF          = 2*pi/kf_mean;     % Chosen spatial scale (Farday wavelength)

%% Dimensionless groups
Reynolds    = xF.^2/(TF*nu);        % Reynolds number
nu0         = 1/Reynolds;           % Inverse Reynolds number
Bo          = sig*TF^2/(rho*xF^3);
G           = g0*TF^2/xF;
M           = drop_mass./(rho*xF^3); % 
cf_air      = 6*pi*drop_radius*mu_air*TF/drop_mass; % Air drag (vector valued for several drops)
c4          = 0.17; % Coefficient of restitution (Molacek)
cf_impact   = c4*sqrt(rho*drop_radius/sig)*TF*g0; % Dissipation during impact

%% Dimensionless depth
h0_deep     = h_top_grid./xF;
d0_deep     = d_deep./xF;
a0_deep     = 0;

d0=  b0/ xF;
%% Dimensionless Faraday wavenumber
kf0_deep     = kf_mean*xF;                   

%% Time dependent part of wave speed, denoted in the notes by \tilde{g}.
% <<< MATT >>> changed the definition of vertical vibration to be
% consistent with my previous research (JFM 2020), for familiarity
% g = @(t) G*(1 + Gam*cos(4*pi*t-theta)); % <<< OLD VALUE >>>
g = @(t) G*(1 - Gam * cos(4*pi*t));

%% Time parameters
dt = dt_desired;
impact_interval = 1;
dt   = impact_interval/ceil(impact_interval/dt);
nsteps_impact = impact_interval/dt;

%% Bottom profile
h = h0_deep;
d = d0_deep.*ones(size(xx));
a = 0.*h_top_grid;

%% Set the wave-numbers and matrix for multiplication in 2D (Notice i is already included)
kx  =  2*pi*1i/Lx*[0:Nx/2-1 0 -Nx/2+1:-1];
ky  =  2*pi*1i/Ly*[0:Ny/2-1 0 -Ny/2+1:-1];
[Kx,Ky] = meshgrid(kx,ky); % <<< MD >>> (more efficient)

% Kx  = zeros(Ny,Nx); Ky = zeros(Ny,Nx);
% 
% for i=1:Ny
%     Kx(i,:) = kx;
% end
% for i=1:Nx
%     Ky(:,i) = ky;
% end

% Three lines below are to make the program faster
K2 = Kx.^2 + Ky.^2;
abs_K = sqrt(-K2);
dissMatrix = exp(2*nu0*dt*(K2));        % Dissipation operator in FS
dissMatrix_half = exp(2*nu0*dt/2*(K2)); % Dissipation operator in FS
shift1 = mod(-[1:Nx]+1,Nx)+1;
shift2 = mod(-[1:Ny]+1,Ny)+1;
KxiKy  = Kx+1i*Ky;
KxmiKy = Kx-1i*Ky;

%% Matrix for second derivatives (do not zero-out highest frequency)  <<< MD >>>
kx_deriv =  2*pi*1i/Lx*[0:(Nx/2-1) (-Nx/2):-1];
ky_deriv =  2*pi*1i/Ly*[0:(Ny/2-1) (-Ny/2):-1];
[Kx_deriv,Ky_deriv] = meshgrid(kx_deriv,ky_deriv);
K2_deriv = Kx_deriv.^2 + Ky_deriv.^2;


%% Stuff for B4

K_vec = linspace(kmin,kmax,Nk)';
dk= K_vec(2)-K_vec(1);
K_vec=K_vec+dk;
K2_vec = K_vec.^2;
K3_vec = K_vec.^3;


h_gridsize=dk;
b4_prefactor=-d0*M*G/(2*pi)*h_gridsize;

%% A5 for H

TD=1/(8*pi^2*nu0);
Me = TD/(1-mem);
beta1 = (8*pi^2*(4*nu0^2+d0*Bo) + d0*G)^2 / (16*nu0*pi^2);

beta_func= @(k) - 1/(Me)-beta1*(k-2*pi).^2;


gammaf_dimensional =  g0 * (Gam/mem);
gamma_dimensional = g0 * Gam;

kC= omega0 / (2 *(b0*g0 + sqrt(b0^2*g0^2 + omega0^2*(4*nu^2+b0*sig/rho))))^(1/2);
C1=2*nu*kC^2/gammaf_dimensional;

A6_numerator = @(k_dim) (2*kC^2*(4*nu^2 + b0*sig/rho) + b0*g0)*(k_dim-kC) + nu *kC*C1*(gamma_dimensional-gammaf_dimensional);
phifunc = @(k) -pi/4;%-1/2 *atan2(1,A6_numerator(k/(2*pi)*kf_mean)/(nu*kC*omega0));
% -pi/4
% -1/2 * acot(2*mem / omega0 *(Gam-Gam/mem))

A5_activation = @(t,s,m) (tanh(t*s + m)+1)/2;

H_A5= @(t,k) (- exp(beta_func(k)*t)) .*cos(2*pi*t +phifunc(k))./(4*pi*sin(phifunc(k)));
H_A13 = @(t,s,k) 2/(4*pi) *exp(beta_func(k).*(t-s))*cos(2*pi*t + phifunc(k))*cos(2*pi*s-phifunc(k));



alpha_func= @(k) k.* sqrt(d0 *(k.^2*Bo + G));
H_A14 = @(t,k) exp(-2*nu0*k.^2.*t).*sin(alpha_func(k).*t)./alpha_func(k);

%% Dissipation operator (including highest frequency) <<< MD >>
% dissipation term over half timestep, i.e. exp(-2 * nu * k^2 * dt/2)
D = exp(nu0 * K2_deriv * dt);

%% Parameters for b1 in k space with integration factor
kappa = 2*K2_deriv*nu0;
% D_product = @(eta,etaprime)deal(eta,  eta.*kappa+etaprime);


%% Dimensionless parameters for drop pressure (not used if useDeltaDrop option is chosen)
useDeltaDrop = 1;
phi0 = zeros(size(xx));


%% Parameters for B1 x space
%Laplacian w/ periodic BCs
ex = ones(Nx,1);
Lapx = spdiags([ex -2*ex ex], [-1, 0, 1], Nx, Nx); Lapx(1,end)=1; Lapx(end,1)=1; Lapx=Lapx/(hx*hx);
Ix = speye(Nx);

ey = ones(Ny,1);
Lapy = spdiags([ey -2*ey ey], [-1, 0, 1], Ny, Ny); Lapy(1,end)=1; Lapy(end,1)=1; Lapy=Lapy/(hy*hy); 
Iy= speye(Ny);

Lap2d = kron(Lapx,Iy) + kron(Ix,Lapy);


delta_sigma=hx/2;
%deltafunc2d = @(x,y) exp(-(x.^2+y.^2)/(2*delta_sigma^2))/(delta_sigma^2*(2*pi));
delta_unnormalized = @(x,y) exp(-(x.^2+y.^2)/(2*delta_sigma^2));
deltafunc2d= @(x,y) delta_unnormalized(x,y) / sum(sum(delta_unnormalized(xx,yy))) / (hx*hy);
%sum(sum(deltafunc2d(xx,yy)))*(hx*hy)

Pfunc = @(x,y) d0*M(1)*G*deltafunc2d(x,y);
%Pmat=@(xp,yp) del2( Pfunc(xx-xp,yy-yp),hx,hy);
Pvec =@(xp,yp) Lap2d*( reshape(Pfunc(xx-xp,yy-yp),[],1) );

%%
varList = who;

%initialiste a structure
p = struct;

%use dynamic fieldnames
for index = 1:numel(varList)
    p.(varList{index}) = eval(varList{index});
end