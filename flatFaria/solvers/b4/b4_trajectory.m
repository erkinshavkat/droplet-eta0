function [x_data,y_data,t_data,eta_data,ui,vi] = b4_trajectory(p)

t = p.theta/(4*pi);
x_data = zeros(p.nimpacts,1); y_data = zeros(p.nimpacts,1); t_data = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.nimpacts*p.nsteps_impact); 

xi=p.xi; yi=p.yi; ui=p.ui; vi=p.vi;

H_vec=zeros(size(p.K_vec));
dH_vec= ones(size(p.K_vec));

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   

    for nn=1:p.nsteps_impact 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        b4_eta_compute=@(x,y) p.b4_prefactor *sum(p.K3_vec.* H_vec.*besselj(0,p.K_vec.* sqrt((x-xi).^2 +(y-yi).^2 )));

        eta_data(:,:,(n-1)*p.nsteps_impact + nn ) = arrayfun(b4_eta_compute,p.xx,p.yy);
    end

end
end