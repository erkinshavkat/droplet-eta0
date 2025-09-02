clear; close all; clc
addpath(genpath(pwd))
gamma=4.4974;
H=0.005;

% gamma=5.6075;
% H=0.001

% Lx=32;
% Nx=1024;
Lx=16;
Nx=256;
eta0=zeros(Nx,Nx);
Nk=Nx*2; Lk=4*pi;
theta=0.5;
mem=0.98;

% num_waveheight=[];
% activate_waveheight=[];
fig = figure('Position', [0, 0, 1500, 900]); 

eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk,theta,mem);
p.xi = 0; p.yi = 0; p.ui= 0; p.vi = 0;
p.nimpacts = 100;


t = p.theta/(4*pi);




% phi = p.phi0; 
% eta = p.eta0; 
% phi_hat = fft2(phi); 
% eta_hat = fft2(eta);
% xi = p.xi; yi = p.yi; ui = p.ui; vi = p.vi;

x_data= zeros(p.nimpacts,1);
y_data= zeros(p.nimpacts,1);


H_num=zeros(length(p.K_vec),p.nimpacts );
dH_num= zeros(length(p.K_vec),p.nimpacts );


% faria_ax = plot(p.y,zeros(p.Nx,1),'Linewidth',2);
% b4_ax = plot(p.y,zeros(p.Nx,1),'LineWidth',2);
% legend('Faria','B4')

trajy=zeros(p.nimpacts,1);%linspace(-2,2,p.nimpacts);%
dtraj=trajy(2)-trajy(1);

% HkF_num=zeros(p.nimpacts*p.nsteps_impact,1);
% Hkf_formula=zeros(p.nimpacts*p.nsteps_impact,1);
% Hkf_activate=zeros(p.nimpacts*p.nsteps_impact,1);
% init_t=t;
% ts=(1:p.nimpacts*p.nsteps_impact)*p.dt+init_t;
% for i=1:p.nimpacts*p.nsteps_impact
%     [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);
%     Hkf_activate(i)=p.H_A14(t-init_t,2*pi)+p.H_A13(t,init_t,2*pi)*p.A5_activation(t-init_t,1/(2*p.nu0*4*pi^2));
%     Hkf_formula(i)=p.H_A14(t-init_t,2*pi)+p.H_A13(t,init_t,2*pi);
%     Hkf_num(i) = interp1(p.K_vec,H_num,2*pi*1.015);
%     t = t + p.dt;
% end

% plot(ts,Hkf_formula,'LineWidth',2);hold on
% plot(ts,Hkf_activate,'LineWidth',2)
% plot(ts,Hkf_num,'LineWidth',2)
% legend({'A13+A14',' A13*tanh(-2\nuk_f^2 t)+A14','Numerical'})
% title(sprintf('H(k_f,t) k_f rescaled, mem=%.2f',mem))
% xlabel('t/T_F')
% ylabel("H(k_F)")
% ylim([-0.3,0.3])
% xlim([0,p.nimpacts])
% figure;

H_A14=zeros(length(p.K_vec),p.nimpacts);
H_A13=zeros(length(p.K_vec),p.nimpacts);
colors = winter(p.nimpacts); % Colormap for H_A14
hold on;

for n=1:p.nimpacts
    xi=0;yi=trajy(n);
    % disp(['Impact number: ' num2str(n)])
    x_data(n) = xi;    y_data(n) = yi;

    % [ui, vi, phi_hat] = drop_impact_matt(xi, yi, ui, vi, phi_hat, eta_hat, p);

    dH_num(:,n) = 1; 
    impact_ts(n)=t;

    for nn=1:p.nsteps_impact 

        % [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        %[b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);

        t= t+p.dt;
    end

    if mod(n,10)==0
        next_x = xi ; next_y=yi+dtraj;
        
        for impact = 1:n
            s = impact_ts(impact);
            elapsed_time = t - s;
            H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
            H_A13(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)).*p.H_A13(t,s,p.K_vec,p.phifunc);

        end
        b4_eta_compute = @(x, y, impact, H) p.b4_prefactor * sum(p.K3_vec .* H(:,impact) .* besselj(0, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        eta_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_eta_compute(x, y, impact, H), 1:n));

        b4_detady_compute = @(x, y, impact, H) ...
                            (y - y_data(impact)).^2./ sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 ).* ...
                            p.b4_prefactor * sum(p.K_vec.^4 .* H(:,impact) .* ...
                            besselj(1, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        detady_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_detady_compute(x, y, impact, H), 1:n));
        
        % eta=real(ifft2(eta_hat));
        % faria_eta= eta(:,p.Nx/2+1);

        eta_num = arrayfun(@(y) eta_centerline(0, y, H_num), p.y);
        eta_formula = arrayfun(@(y) eta_centerline(0, y, H_A14+H_A13), p.y);

        plot(p.y, eta_num, 'Color', colors(n,:), 'LineWidth', 2, 'DisplayName', sprintf('num H, n=%d', n));
        plot(p.y, eta_formula, 'Color', colors(n,:), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', sprintf('A14+A13 var \phi, n=%d', n));
    end
end
xlim([-4 4])
xlabel('y');
ylabel('\eta');
title(['\eta(0,y) at impact, \theta=',num2str(p.theta/pi),'\pi, k_C=kf_mean']);
legend('show');
saveas(fig, sprintf('vis/eta A13 varphase theta%.2fpi kC=kf_mean.png',p.theta/pi));