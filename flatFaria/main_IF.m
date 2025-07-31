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


toc
