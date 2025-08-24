function [H_A5, H_A14] = optimize_activation(p)


H_vec=zeros(length(p.K_vec),1);
dH_vec= ones(length(p.K_vec),1);

N=p.nimpacts*p.nsteps_impact;

H_num_data=zeros(p.Nk,N);
t_vec=(0:N-1)*p.dt;
for n=1:N
    t=t_vec(n);
    H_num_data(:,n) = H_vec;
    [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
end

activ=p.A5_activation;

H_formula = @(t, k,tau) activ(t,tau) * p.H_A13(t,0, k) + p.H_A14(t, k);

[t_domain ,k_domain] = meshgrid(t_vec, p.K_vec);

compute_H_formula = @(tau) arrayfun(@(t, k) H_formula(t, k, tau), t_domain, k_domain);
error = @(tau) norm(H_num_data - compute_H_formula(tau), 'fro');
mins=fminsearch(@(tau) error(tau), 1);
mins
end





