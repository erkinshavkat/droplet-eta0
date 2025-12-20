clear; close all; clc
addpath(genpath(pwd))
gamma=4.4974;

Lx=16;
Nx=256;
eta0=zeros(Nx,Nx);
Nk=Nx; Lk=6*pi;
eta0=zeros(Nx);
mem=1;
theta=0;

H_vals=[0.0005,0.001,0.002,0.003,0.005,0.01];
omega0_vals=[80,160,400,800];
nH = length(H_vals);
nOmega = length(omega0_vals);

fig = figure('Position', [0, 0, nOmega*700, nH*500]);

for h_idx=1:nH
    for omega_idx=1:nOmega
        H = H_vals(h_idx);
        omega0 = omega0_vals(omega_idx);
        p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk,theta,mem,omega0);

        x_domain=linspace(-5,5,Nx/2);
        
        subplot(nH, nOmega, (h_idx-1)*nOmega + omega_idx);
        ax=plot(x_domain,0*x_domain,"LineWidth",4,"Color",'black',"DisplayName",'k tanh(kh)(Reference)'); hold on
        ax2=plot(x_domain,0*x_domain,"LineWidth",2,"DisplayName",'hk^2 (Shallow)'); hold on
        ax3=plot(x_domain,0*x_domain,"LineWidth",2,"DisplayName",'bk^2 (Faria)'); hold on
        ax4=plot(x_domain,0*x_domain,"LineWidth",2,"DisplayName",'k (Deep/Milewski)'); hold on

        DtNs=[p.K_vec.*tanh(p.K_vec.*p.h0) p.K2_vec.*p.h0  p.K2_vec.*p.d0 p.K_vec]; 
        axes={ax,ax2,ax3,ax4};
        
        for t=100.99
            impacts=floor(t);
            range = 0;
            for j=1:4
                DtN=DtNs(:,j);
                ax=axes{j};
                HL=zeros(Nk,1);
                for i=0:impacts
                    HL = HL+HL_varDtN(t-i,p.K_vec,p,DtN);
                end
                etaL=zeros(Nx/2,1);
                for i=1:length(x_domain)
                    x = x_domain(i);
                    etaL(i) = eta_varDtN(x,0,0,0,HL,DtN,p);
                end
                ax.YData = etaL;
                if j==1
                    range = max(abs(etaL))*1.5;
                end
            end
            ylim([-range,range])
            xlim([-5,5])
            legend()
            title(sprintf('h=%.2fmm, \\omega=%dHz',H*1000,p.omega0/(2*pi)))
        end
    end
end

saveas(gcf,'vis/DtN_comp/all_comparisons.png');