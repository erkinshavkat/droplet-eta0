clear 
clc 
addpath('solvers')

% The threshold test should be run with the SAME Nx and Lx as the
% simulations of interest. Also make sure that the geometry is the same
Nx=168;
Lx=18;
H=0.005;
eta0=zeros(Nx);
% Gam = 3.8:0.01:3.91;
Gam = [3.6,3.65,3.7];
gamma_n = 1;
%%
for m=1:gamma_n
    for n=1:length(Gam)
        p = setup_IF_matt(Gam(n),H,eta0,Nx,Lx,n);
        % Number of impacts (time to wait to check the threshold)
        p.nimpacts = 300; % 1000
        [is_stable(n),eta_max(n),~,sqrt_energy(n,:),t_vec(n,:)] = faraday_threshold(p,[0 0]); 
    end
       save([pwd, '/results/' sprintf('%.5f',m) '_faraday.mat'],'Gam','eta_max','sqrt_energy','t_vec');
end

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
