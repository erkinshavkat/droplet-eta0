clear; close all; clc
addpath(genpath(pwd))
gamma=4.4974;
H=0.005;

% gamma=5.6075;
% H=0.001

% Lx=32;
% Nx=1024;
Lx=16
Nx=256
eta0=zeros(Nx,Nx);
Nk=Nx; Lk=12*pi;
theta=1.2;


num_waveheight=[];
activate_waveheight=[];
fig = figure('Position', [0, 0, 1500, 900]); 
thetas=linspace(0,2,50);


for theta=thetas

eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk,theta);
p.xi = 0; p.yi = 0; p.ui= 0; p.vi = 0;
p.nimpacts = 30;


t = p.theta/(4*pi);




phi = p.phi0; 
eta = p.eta0; 
phi_hat = fft2(phi); 
eta_hat = fft2(eta);
xi = p.xi; yi = p.yi; ui = p.ui; vi = p.vi;

x_data= zeros(p.nimpacts,1);
y_data= zeros(p.nimpacts,1);


H_num=zeros(length(p.K_vec),p.nimpacts);
dH_num= zeros(length(p.K_vec),p.nimpacts);




t_vec=[];

trajy=zeros(p.nimpacts,1);%linspace(-2,2,p.nimpacts);%
dtraj=trajy(2)-trajy(1);

num_taus=15;
taus=linspace(0,1.5,num_taus);
H_A14=zeros(length(p.K_vec),p.nimpacts);
H_A13_activate=zeros(length(p.K_vec),p.nimpacts,num_taus);



for n=1:p.nimpacts
    xi=0;yi=trajy(n);
    disp(['Impact number: ' num2str(n)])
    x_data(n) = xi;    y_data(n) = yi;

    [ui, vi, phi_hat] = drop_impact_matt(xi, yi, ui, vi, phi_hat, eta_hat, p);

    dH_num(:,n) = 1; 
    impact_ts(n)=t;

    for nn=1:p.nsteps_impact 

        %[phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        %[b1k.eta_hat, b1k.etaprime_hat] = b1k_evolve_wave_rkstep(b1k.eta_hat,b1k.etaprime_hat, t + (nn -1)*p.dt, p); 
        [H_num, dH_num] = H_eq_rkstep(H_num,dH_num, t, p);

        t= t+p.dt;


        % if n>28
        %     next_x = xi ; next_y=yi+dtraj;
        %     for impact = 1:n
        %         s = impact_ts(impact);
        %         elapsed_time = t - s;
        %         H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
        %         for nnn=1:num_taus
        %             H_A13_activate(:,impact,nnn) = p.A5_activation(elapsed_time,taus(nnn)).*p.H_A13(t,s,p.K_vec);
        %         end
        %     % disp([next_x,next_y])
        %     % disp(sqrt((next_x - x_data(impact)).^2 + (next_y - y_data(impact)).^2 ))
        %     end
        %     b4_eta_compute = @(x, y, impact, H) p.b4_prefactor * sum(p.K3_vec .* H(:,impact) .* besselj(0, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        %     eta_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_eta_compute(x, y, impact, H), 1:n));

        %     b4_detady_compute = @(x, y, impact, H) ...
        %                         (y - y_data(impact)).^2./ sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 ).* ...
        %                         p.b4_prefactor * sum(p.K_vec.^4 .* H(:,impact) .* ...
        %                         besselj(1, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        %     detady_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_detady_compute(x, y, impact, H), 1:n));
            
            
        %     eta00_num = eta_centerline(next_x,next_y,H_num);
        %     eta00_formula = arrayfun(@(n) eta_centerline(next_x,next_y,H_A13_activate(:,:,n)) , 1:num_taus);
        %     num_waveheight=([num_waveheight,eta00_num]);
        %     activate_waveheight=([activate_waveheight;eta00_formula]);
        %     t_vec=[t_vec,t];
        % end

    end

end
next_x = xi ; next_y=yi+dtraj;
for impact = 1:n
    s = impact_ts(impact);
    elapsed_time = t - s;
    H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
    for nnn=1:num_taus
        H_A13_activate(:,impact,nnn) = p.A5_activation(elapsed_time,taus(nnn)).*p.H_A13(t,s,p.K_vec);
    end
end
b4_eta_compute = @(x, y, impact, H) p.b4_prefactor * ...
    sum(p.K3_vec .* H(:,impact) .*...
    besselj(0, p.K_vec .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
eta_centerline = @(x,y,H) sum(arrayfun(@(impact) b4_eta_compute(x, y, impact, H), 1:n));
eta00_num = eta_centerline(next_x,next_y,H_num);
eta00_formula = arrayfun(@(n) eta_centerline(next_x,next_y,H_A13_activate(:,:,n)) , 1:num_taus);
num_waveheight=([num_waveheight,eta00_num]);
activate_waveheight=([activate_waveheight;eta00_formula]);

end
disp(num_waveheight)
disp(activate_waveheight)
plot(thetas, num_waveheight,"LineWidth",3,"Color",'red');hold on
legend_labels = {'num'};
title('\eta(0,0,t), bouncer, activation tanh(t/\tau)');

% Create a color gradient for the tau curves
cmap = winter(num_taus);

for n=1:num_taus
    plot(thetas, activate_waveheight(:,n), 'LineWidth', 3, 'Color', cmap(n,:));
    legend_labels{end+1} = ['\tau=',sprintf('%.2f', taus(n)),'T_F'];
end

legend(legend_labels)
xlabel('\theta_I/\pi')
ylabel('\eta')

% disp(eta00_formula)
% plot(t_vec, num_waveheight,"LineWidth",3,"Color",'red');hold on
% legend_labels = {'num'};
% title(['\eta(0,0,t), \theta_I=',sprintf('%.2f', p.theta/pi),'\pi, with activation tanh(t/\tau)']);

% % Create a color gradient for the tau curves
% cmap = winter(num_taus);

% for n=1:num_taus
%     plot(t_vec, activate_waveheight(:,n), 'LineWidth', 3, 'Color', cmap(n,:));
%     legend_labels{end+1} = ['\tau=',sprintf('%.2f', taus(n)),'T_F'];
% end

% for i=28:p.nimpacts
%     xline( p.theta/(4*pi) + i ,"LineStyle",':','LineWidth',2)
% end
% legend(legend_labels)
% xlabel('t/T_F')
% ylabel('\eta')

% ylim([-0.02 0.02])