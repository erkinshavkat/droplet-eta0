clear
close all
mkdir('results') %Create output directory
addpath(genpath(pwd))
%% Initial Condition Setup


tic
gamma=3.5588;
H=0.005;
Lx=16;
Nx=128;

eta0=zeros(Nx,Nx);
Nk=Nx; Lk=Lx;

hole_IF(0,0,0,0,gamma,H,Nx,Lx,Nk*2,Lk*2);

% p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk);
% rho=p.rho; g=p.g0; sig=p.sig; b0=p.b0; rb0=sqrt(p.b0); omega0=p.omega0; nu=p.nu;
% gamma_formula = 4*nu*omega0/ b0 /g
% kF_formula = rho*omega0^2/(2*rb0*(rho*g*rb0+sqrt(sig*rho*omega0^2+b0*rho^2*g^2)))
% kF_sim=p.kf_mean^2
toc
