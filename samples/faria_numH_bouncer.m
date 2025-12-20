%% example comparing Faria solver and wavefield reconstructed with numerical H
clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=128;
Nk=Nx; Lk=12*pi;
mem=0.98; theta=1.3;
eta0=zeros(Nx);
omega=80;
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);


nimpacts=100;


t = p.theta/(4*pi);

%setup H code
H_vec=zeros(Nk,nimpacts);
dH_vec=zeros(Nk,nimpacts);

%setup Faria code
phi = p.phi0; eta = p.eta0; 
xi=0; yi=0; ui=0; vi=0; 
phi_hat = fft2(phi);eta_hat = fft2(eta);


%plotting
faria_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","Faria"); hold on
H_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","H"); hold on
for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])
    %Faria impact by computing delta jump for phi_hat
    [~, ~, phi_hat] = drop_impact_matt(xi,yi, ui, vi, phi_hat, eta_hat, p);

    %H impact by setting dH=1
    dH_vec(:,n)=1;
    
    for nn=1:p.nsteps_impact 
        % loop for time integration between impacts
        
        %integrate Faria and H
        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec, dH_vec, t, p);

        %compute Faria eta with ift
        faria_eta=real(ifft2(eta_hat));

        %evaluate eta from H. b4_eta can only evaluate one point at once. To evaluate along (x,0), need loop
        H_eta=zeros(Nx,1);
        for index=1:Nx
            H_eta(index)=b4_eta(p.x(index),0,0,0,H_vec,p);
        end

        faria_ax.YData=faria_eta(Nx/2+1,:);
        H_ax.YData=H_eta;
        pause(0.01)
        t=t+p.dt;
    end



end