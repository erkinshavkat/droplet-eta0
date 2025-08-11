clear; close all; clc
addpath(genpath(pwd))
gamma=3.5588;
H=0.005;
Lx=16;
Nx=128;

Nk=Nx; Lk=Lx;
eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk);
p.xi = 0; p.yi = 0; p.ui= 0; p.vi = 0;
p.nimpacts = 10;


t1= 0;
H1=zeros(length(p.K_vec),1);
dH1= ones(length(p.K_vec),1);

t2=1.3/4; %theta_I/omega_0
H2=zeros(length(p.K_vec),1);
dH2= ones(length(p.K_vec),1);


%faria_ax=plot(p.y,zeros(p.Nx,1),'LineWidth',2);hold on;
H1_ax=plot(p.K_vec/(2*pi),H1,'LineWidth',2);hold on;
H2_ax=plot(p.K_vec/(2*pi),H2,':',"LineWidth",2);
xlabel('k/kF'); ylabel('H');
legend('H1','H2')


for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])

    for nn=1:p.nsteps_impact 

        [H1, dH1] = H_eq_rkstep(H1,dH1, t1, p);
        [H2, dH2] = H_eq_rkstep(H2,dH2, t2, p);


        if mod(nn,2)==0
            H1_ax.YData=H1;
            H2_ax.YData=H2;


            ylim([-0.5 0.5])


            title(sprintf('mem=%.2f A5+A14, t=%f Tf', p.mem, t1));
            pause(1/12); 
        end
        t1= t1+p.dt;
        t2= t2+p.dt;
    end


end