function [x_data,y_data,eta_data,eta_intermediate] = b4_trajectory(p)

t = p.theta/(4*pi);
x_data = zeros(p.nimpacts,1); y_data = zeros(p.nimpacts,1); t_data = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 
eta_intermediate= zeros(p.Nx,p.Ny,p.nsteps_impact); 


xi=p.xi; yi=p.yi; ui=p.ui; vi=p.vi;

H_vec=zeros(length(p.K_vec),p.nimpacts);
dH_vec= zeros(length(p.K_vec),p.nimpacts);

for n=1:p.nimpacts

    disp(['Impact number: ' num2str(n)])

    x_data(n) = xi;    y_data(n) = yi;    t_data(n) = t;
    dH_vec(:,n) = 1;

    for nn=1:p.nsteps_impact
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);
        if n>p.nimpacts-1
            eta_intermediate(:,:,n) = arrayfun(@(x,y) b4_multiple_impact(x,y,x_data,y_data,H_vec,p),p.xx,p.yy);

        end
    end

    eta_data(:,:,n) = arrayfun(@(x,y) b4_multiple_impact(x,y,x_data,y_data,H_vec,p),p.xx,p.yy);

end
end