%% example comparing etaF computed with HF variable phi, HF constant phi, and the approximate etaF formula
clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=256;
Nk=Nx; Lk=12*pi;
mem=0.98; theta=1.3;
eta0=zeros(Nx);
omega=80;
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);


nimpacts=100;


t = p.theta/(4*pi);

varphi_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","HF variable \phi"); hold on
constphi_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","HF constant \phi"); hold on
etaFapprox_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","etaF approximate formula"); hold on
legend()
xlim([-3 3])
ylim([-0.02 0.02])

%preallocate HFs
Hvarphi=zeros(Nk,nimpacts);
Hconstphi=zeros(Nk,nimpacts);


impact_times=[];
for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])

    %recording impact times
    impact_times=[impact_times t];

    steps=10;
    for nn=1:steps
        % loop for time between impacts
        % Here the time step can be whatever since we are no longer integrating
        % just make sure to adjust dt 
  

        etaF_formula=zeros(1,Nx);
        for nnn=1:length(impact_times)
            s=impact_times(nnn);
            %computing HF from each impact, storing as columns
            Hconstphi(:,nnn)=HF_formula(t,s,p.K_vec,p);
            Hvarphi(:,nnn)=HF_formula(t,s,p.K_vec,p,'cotan');
            
            etaF_formula=etaF_formula+etaF_approx(p.x,0,t,s,p);
        end

        eta_varphi=zeros(Nx,1);
        eta_constphi=zeros(Nx,1);

        for index=1:Nx
            eta_varphi(index)=b4_eta(p.x(index),0,0,0,Hvarphi,p);
            eta_constphi(index)=b4_eta(p.x(index),0,0,0,Hconstphi,p);
        end

        varphi_ax.YData=eta_varphi;
        constphi_ax.YData=eta_constphi;
        etaFapprox_ax.YData=etaF_formula;
        pause(0.01)
        t=t+1/steps;
    end

end