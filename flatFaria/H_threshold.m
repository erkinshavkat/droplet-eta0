clear 
clc 
addpath('solvers')

% The threshold test should be run with the SAME Nx and Lx as the
% simulations of interest. Also make sure that the geometry is the same
Nx=1024;
Lx=32;
H=0.005;
Nk=Nx;
eta0=zeros(Nx);
% Gam = 3.8:0.01:3.91;
Gams = linspace(4.5,4.6,5);

% Create figure for animation of (p.K_vec, H_vec)


growth_rates= zeros(length(Gams), 1);
for i = 1:length(Gams)
    i
    gamma = Gams(i);
    p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,2*pi*0.9,2*pi*1.1,1.1);
    p.nimpacts = 300;
    H_vec = zeros(length(p.K_vec),1);
    dH_vec = ones(length(p.K_vec),1);
    t = 0;
    peak_values = zeros(p.nimpacts,1);

    % figure_anim = figure;
    % h_anim = plot(p.K_vec, H_vec, 'LineWidth', 2);
    % xlabel('k'); ylabel('H'); title('H vs k (Animation)');
    % grid on;
    % drawnow;
    for n = 1:p.nimpacts
        max_H = 0;
        for nn = 1:p.nsteps_impact

            [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
            max_H = max([max(abs(H_vec)), max_H]);
            t= t + p.dt;
            % figure(figure_anim);
            % set(h_anim, 'YData', H_vec);
            % title(sprintf('H vs k, Impact %d', n));
            % drawnow;
        end
        peak_values(n) = max_H;

        % Update animation figure



    end

    % Plot peak_values for this gamma in the peak figure

    % semilogy(1:p.nimpacts, peak_values); hold on;
    semilogy(peak_values, 'o-'); hold on;
    m=30;
    logvalue=log(peak_values);
    rate1= (logvalue(end) - logvalue(end-1));
    rate2= (logvalue(end) - logvalue(end-m))/m;
    if abs(1-rate2/rate1) < 0.1
        growth_rates(i) = rate1;
    else
        growth_rates(i) = NaN;
    end


end
A=[ones(length(Gams),1) Gams.'];
b=growth_rates;
x = (A.'*A)\(A.' * b);
a = x(1); b = x(2);
disp(-a/b)
plot(Gams, growth_rates, 'o-');


%% Estimate Faraday threshold
% Let nu be the growth/decay rate of the wave field, where nu > 0 for
% growth (above the Faraday threshold) and nu < 0 for decay (below the
% Faraday threshold). Our goal is to use linear interpolation to estimate
% the value of Gamma, denoted Gamma_F, for which nu = 0 (the Faraday
% threshold). This formulation assumes that the simulation is run for long
% enough for there to be exponential growth/decay in nu over time, which
% may be hard to achieve near the Faraday threshold. We thus exclude any
% values of Gamma in our estimation of the Faraday threshold that are too
% close to the Faraday threhsold (perhaps counter-intuitively) as these
% values are too unreliable -- it is very hard to tell whether there will
% be growth or decay (and at what rates) when there is too little data.
% Running longer simulations may make these data points more useful.

% Step 1. Estimate the value of nu for each value of Gamma, and exclude any
% values of Gamma for which nu does not appear to have approached
% exponential growth/decay in the simulation window.

% Tolerance for checking that the growth is sufficiently linear near the
% end of the simulation. We compute the slope of log_energy (define below)
% using the final Faraday period (estimate 1) and the final N Faraday
% periods (estimate 2). If these two values are sufficiently close (within 
% TOL), then we use estimate 1. If these two value are far apart (the
% difference is >= TOL) then we discard this value of Gamma.

% Determining the correct tolerance (and value of N, below) may take some
% fine-tuning. Another option is to just look at the data and remove any
% values of Gamma (Step 2) in which we do not see a straight line in the
% following plot for large t/T_F.
return 
TOL = 1e-1; 
N = 30; % final N Faraday periods used for estimate 2. I suggest N = 30;

nu_estimate_vec = zeros(1, length(Gam));

figure; 
hold on
for n = 1:length(Gam)
    % Take logaritham of growth rate over time. If nu has exponential
    % growth/decay then log_energy has linear growth/decay.
    log_energy = log(sqrt_energy(n,:)); 
    
    % Plot data to ensure that there is long-time linear behaviour in
    % log_energy.
    plot(t_vec(n,:), log_energy, 'LineWidth', 2); 
    
    % Estimate growth/decay rate at final time using a backward difference
    estimate_1 = log_energy(end) - log_energy(end-1);
    estimate_2 = (log_energy(end) - log_energy(end-N))/N;
    abs(1-estimate_2/estimate_1)
    if abs(1-estimate_2/estimate_1) < TOL
        nu_estimate_vec(n) = estimate_1;
    else
        nu_estimate_vec(n) = NaN;
    end
end

set(gca, 'FontSize', 20, 'FontName', 'Times New Roman')
xlabel('$$t/T_F$$', 'Interpreter', 'Latex', 'FontSize', 22); 
ylabel('$$\log(E)$$', 'Interpreter', 'Latex', 'FontSize', 22);
legend_cell = cell(1, length(Gam));
for ii = 1:length(Gam); legend_cell{ii} = ['$$\Gamma = $$ ' num2str(Gam(ii))]; end
legend(legend_cell, 'Interpreter', 'Latex', ...
    'FontSize', 22, 'Location', 'Best')
grid on
box on

% Step 2. Perform linear interpolation of the vector nu_estimate_vec over
% the data Gam, excluding any data points for which nu_estimate_vec = NaN.
ind = ~isnan(nu_estimate_vec);
nu_estimate_sub = nu_estimate_vec(ind);
Gam_sub = Gam(ind);

% create Vandermonde matrix for linear least squares. We want to minimise
% ||A*x - b|| in the 2-norm, where x is the vector of parameters for the
% linear fit.
A = [ones(length(Gam_sub),1) Gam_sub.'];
b = nu_estimate_sub.';

% compute the parameters by solving the normal equations, 
% (A.'* A) * x = A.' * b. We could also use the QR-decomposition or SVD.
x = (A.'*A)\(A.' * b);

% Decompose parameters, with nu(Gamma) = a + b * Gamma
a = x(1); b = x(2);

% Find Faraday threshold, so that nu(Gamma_F) = 0
Gamma_F_est = -a/b;

% legend off
% grid off

axis square
% exportgraphics(gcf,'C:\Users\abeja\OneDrive - University of North Carolina at Chapel Hill\Anderson Localization\Figures\Supplementary\SFigure7\Materials\logE.pdf','ContentType','vector');

%%
figure
hold on
fplot(@(Gam) a + b * Gam, [min(Gam) max(Gam)], 'b', 'LineWidth', 2)
plot(Gam, nu_estimate_vec, 'ko', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
plot(Gamma_F_est, 0, 'rd', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
grid on
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
xlabel('$$\Gamma = \gamma/g$$', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel('Growth rate', 'Interpreter', 'Latex', 'FontSize', 22)
title(['$$\Gamma_F = $$ ' num2str(Gamma_F_est)], ...
    'Interpreter', 'Latex', 'FontSize', 22)

axis square
% exportgraphics(gcf,'C:\Users\abeja\OneDrive - University of North Carolina at Chapel Hill\Anderson Localization\Figures\Supplementary\SFigure7\Materials\fit.pdf','ContentType','vector');
