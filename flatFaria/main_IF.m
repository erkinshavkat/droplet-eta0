clear
close all
mkdir('results') %Create output directory
addpath(genpath(pwd))
%% Initial Condition Setup


tic
gamma=3.5588;
H=0.005;
Nx=96;
Lx=16;
eta0=zeros(Nx,Nx);

hole_IF(0,0,0.05,0,gamma,H,eta0,Nx,Lx,1);

toc