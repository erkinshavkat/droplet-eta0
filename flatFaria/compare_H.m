function [eta_b4, eta_faria] = compare_H(p)
t = p.theta/(4*pi);

H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

fig = figure('Position', [0, 0, 1500, 900]); 

K_scales = [1,1.015,1.05];
num_scales = numel(K_scales);
H_formula_axes = gobjects(num_scales,1);
legend_labels = cell(num_scales,1);

for kidx = 1:num_scales
    K_vec_scaled = p.K_vec * K_scales(kidx);
    H_formula_axes(kidx) = plot(K_vec_scaled/(2*pi), zeros(p.Nk,1), 'LineWidth', 2); hold on;
    legend_labels{kidx} = sprintf('A13+A14, K_{scale}=%.3f', K_scales(kidx));
end
H_activation_ax = plot(p.K_vec/(2*pi)*1.015, zeros(p.Nk,1),'LineWidth', 2);
legend_labels{end+1} = 'A13*tanh(t/\tau)+A14, K_{scale}=1.015';

H_num_ax = plot(p.K_vec/(2*pi), zeros(p.Nk,1), 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'numerical';



xline(1,'k:','LineWidth',2);

xlabel('k/k_F'); ylabel('H');
legend(legend_labels);

L2errs=zeros(p.nimpacts*p.nsteps_impact,1);

% v = VideoWriter(sprintf('vis/A13 tanh vs no tanh theta %.2fpi.avi',p.theta/pi),'Motion JPEG AVI');
% v.FrameRate = 24;
% open(v);

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])

    dH_vec(:,n)=1;
    impact_ts(n)=t;
    for nn=1:p.nsteps_impact 

        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        t= t+p.dt;

        if n>p.nimpacts-3
            for impact = 1:n
                s = impact_ts(impact);
                elapsed_time = t - s;
                H_A14(:,impact) = p.H_A14(elapsed_time,p.K_vec);
                H_A13(:,impact) = p.H_A13(t,s,p.K_vec);
                H_A13_activation(:,impact) = p.A5_activation(elapsed_time,0.7775) *p.H_A13(t,s,p.K_vec);
            end
            sum_H_num = sum(H_vec,2);
            sum_H_formula = sum(H_A13+H_A14,2);
            sum_H_activate = sum(H_A13_activation+H_A14,2);

            % Update all formula axes with sum_H_formula
            for kidx = 1:num_scales
                H_formula_axes(kidx).YData = sum_H_formula;
            end
            H_num_ax.YData = sum_H_num;
            H_activation_ax.YData=sum_H_activate;

            % frame = getframe(gcf);
            % writeVideo(v, frame);

            L2err = norm(sum_H_num - sum_H_formula);
            L2errs((n-1)*p.nsteps_impact + nn) = L2err;

            ylim([-4 4])

            title(sprintf('mem=%.2f, t=%f Tf, theta=%.2fpi', p.mem, t,p.theta/pi));

            pause(1/24); 
        end
    end
end

% close(v);
end





