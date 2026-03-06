%% example comparing etaL computed with exact HL, exact HLinf, and the approximate etaL formula
clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;
fig = figure('Position', [100, 100, 1200, 800]); 

Lx=16; Nx=256;
Nk=Nx; Lk=12*pi;
mem=0; theta=1.3;
eta0=zeros(Nx);

omega=80;

p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);


nimpacts=10;


t0 = p.theta/(4*pi);
t=t0;
eta_faria_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","Exact \eta_L discrete"); hold on
eta_discrete_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","discrete formula"); hold on
eta_cont_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","cont formula"); hold on
eta_nodelta_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","gaussian impact R=4, delta->cos"); hold on

ylim([-0.01 0.01])
xlim([-6 6])
legend('FontSize',18)




impact_times=[];

phi = p.phi0; eta = p.eta0; 
xi=0; yi=0; ui=0; vi=0; 
phi_hat = fft2(phi);eta_hat = fft2(eta);

%maintain speed of 5mm/s 
speed=5;
motionrange = (speed/1000)/p.xF *p.TF * nimpacts

% Create video writer
videoFile = sprintf('vis/cont_widegaussianimpact_approxdelta_walk%.dmm_omega%.d.avi',speed,omega);
v = VideoWriter(videoFile);
v.FrameRate = 18;
open(v);


xtraj=linspace(-motionrange/2,motionrange/2,nimpacts); 
ytraj=zeros(nimpacts,1);

xtraj_fine=linspace(-motionrange/2,motionrange/2,nimpacts*p.nsteps_impact);
ytraj_fine=zeros(nimpacts*p.nsteps_impact,1);

R=4*p.drop_radius / p.xF;
for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])
    tic
    %recording impact times
    [~, ~, phi_hat] = drop_impact_matt(xtraj(n),0, ui, vi, phi_hat, eta_hat, p);

    impact_times=[impact_times t];
    for nn=1:p.nsteps_impact 
        % loop for time between impacts
        % Here the time step can be whatever since we are no longer integrating
        % just make sure to adjust dt 
        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 

        if (mod(nn,8)==0) %plotting, only picking every 8 steps
            eta_discrete=zeros(1,Nx);eta_strobe=zeros(1,Nx);

            for nnn=1:length(impact_times)
                s=impact_times(nnn);


                eta_discrete=eta_discrete+etaL_approx(p.x-xtraj(nnn),0,t-s,p);
            end
            eta_cont=zeros(1,Nx);eta_nodelta=zeros(1,Nx);
            for nnn=1:round(((t)-t0)/p.dt)
                s=t0+(nnn-1)*p.dt;
                eta_nodelta=eta_nodelta+ p.dt.*etaL_gaussian(p.x-xtraj_fine(nnn),0,t-s,R,p).*(1+2*cos(2*pi*(t-s)));
                eta_cont=eta_cont+ p.dt.*etaL_approx(p.x-xtraj_fine(nnn),0,t-s,p);
            end

            eta0_exact=real(ifft2(eta_hat));
            eta_faria_ax.YData=gather(eta0_exact(Nx/2+1,:));
            eta_discrete_ax.YData=gather(eta_discrete);
            eta_nodelta_ax.YData=gather(eta_nodelta);
            eta_cont_ax.YData=gather(eta_cont);
            title(['Impact: ' num2str(t-t0)])
            drawnow;
            writeVideo(v, getframe(gcf));
        end
        t=t+p.dt;
    end
    %plotting frame freeze before impact
    for nn=1:p.nsteps_impact/10
    writeVideo(v, getframe(gcf));
    end

    toc
end

close(v);