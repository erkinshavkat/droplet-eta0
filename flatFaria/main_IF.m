clear
close all
mkdir('results') %Create output directory
addpath(genpath(pwd))
%% Initial Condition Setup


tic
gamma=3.5588;
H=0.005;
Nx=128;
Lx=16;
eta0=zeros(Nx,Nx);
for Nk=[1,2,4]
hole_IF(0,0,0,0,gamma,H,eta0,Nx,Lx,Nk);
end

for Nx=[256,384,512]
    hole_IF(0,0,0,0,gamma,H,eta0,Nx,Lx,Nk);
end
toc