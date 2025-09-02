function [eta_b4, eta_faria] = compare_H(p)
t = p.theta/(4*pi);

H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1500, 900]); 


H_num_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), 'LineWidth', 2);hold on;
H_varphi_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), '--', 'LineWidth', 2);
H_fixphi_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), '--', 'LineWidth', 2);
H_lin_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), ':', 'LineWidth', 2); % New axis for H_A13_lin
ylim([-1 1])
xlim([0 2])
xline(1,'k:','LineWidth',2);
yline(0,'k:','LineWidth',2);

xlabel('k/k_F'); ylabel('H');
legend('numerical','A13 var \phi', 'A13 fixed \phi=\pi/4', 'A13 \phi linear approx')

v = VideoWriter(sprintf('vis/H linear phi A13 bouncer.avi',p.theta/pi),'Motion JPEG AVI');
v.FrameRate = 24;
open(v);

for n=1:p.nimpacts
    disp(['Impact number: ' num2str(n)])
    dH_vec(:,n)=1;
    impact_ts(n)=t;
    for nn=1:p.nsteps_impact
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        t= t+p.dt;

        if n>p.nimpacts-5
            for impact = 1:n
                s = impact_ts(impact);
                elapsed_time = t - s;
                H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
                H_A13_varphi(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) * p.H_A13(t, s, p.K_vec, p.phifunc);
                H_A13_fixphi(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) * p.H_A13(t, s, p.K_vec, @(k) -pi/4);
                H_A13_lin(:,impact) = p.A5_activation(elapsed_time,1/(2*p.nu0*4*pi^2)) * p.H_A13(t, s, p.K_vec, p.philin); % New
            end
            sum_H_num = sum(H_vec,2);
            sum_H_varphi = sum(H_A14+ H_A13_varphi,2);
            sum_H_fixphi = sum(H_A14+ H_A13_fixphi,2);
            sum_H_lin = sum(H_A14+ H_A13_lin,2); % New

            H_num_ax.YData = sum_H_num;
            H_varphi_ax.YData = sum_H_varphi;
            H_fixphi_ax.YData = sum_H_fixphi;
            H_lin_ax.YData = sum_H_lin; % New

            frame = getframe(gcf);
            writeVideo(v, frame);

            title(sprintf('mem=%.2f, t=%f Tf, theta=%.2fpi', p.mem, t, p.theta/pi));
            pause(1/24); 
        end
    end
end

close(v);
end





