clear

mkdir('results') %Create output directory
addpath(genpath(pwd))
%% Initial Condition Setup
n = 12; %number of simulations %%%% MAIN THING TO EDIT


tic
gamma=3.5588;
H=0.005;
Nx=128;
Lx=16;
eta0=zeros(Nx,Nx);

hole_IF(0,0,0.05,0,gamma,H,eta0,Nx,Lx);

toc