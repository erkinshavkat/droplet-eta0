%% example comparing etaL computed with exact HL, exact HLinf, and the approximate etaL formula
clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=128;
Nk=Nx; Lk=12*pi;
mem=0; theta=1.3;
eta0=zeros(Nx);
omega=80;
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);


nimpacts=100;


t = p.theta/(4*pi);

HL_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","HL"); hold on
HLinf_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","HLinf"); hold on
etaLapprox_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","etaL approximate formula"); hold on
ylim([-0.002 0.002])
legend()
%preallocate HL
HLfunc=zeros(Nk,nimpacts);

impact_times=[];
for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])

    %recording impact times
    impact_times=[impact_times t];
    for nn=1:p.nsteps_impact 
        % loop for time between impacts
        % Here the time step can be whatever since we are no longer integrating
        % just make sure to adjust dt 

        pause(0.01)
        t=t+p.dt;
    end
        etaL_formula=zeros(1,Nx);
        for nnn=1:length(impact_times)
            s=impact_times(nnn);
            %computing HL from each impact, storing as columns
            HLfunc(:,nnn)=HL_formula(t-s,p.K_vec,p);
            
            %computing etaL from each impact and superimpose
            etaL_formula=etaL_formula+etaL_approx(p.x,0,t-s,p);
        end
        %compute HLinf only once since tau denotes time since most recent impact
        HLinf=Hinftau(t-impact_times(end),p.K_vec,p);

        %for HL and HLinf, evaluating the wavefield along (x,0)
        eta_HL=zeros(Nx,1);
        eta_HLinf=zeros(Nx,1);

        for index=1:Nx
            eta_HL(index)=b4_eta(p.x(index),0,0,0,HLfunc,p);
            eta_HLinf(index)=b4_eta(p.x(index),0,0,0,HLinf,p);
        end

        HL_ax.YData=eta_HL;
        HLinf_ax.YData=eta_HLinf;
        etaLapprox_ax.YData=etaL_formula;


end