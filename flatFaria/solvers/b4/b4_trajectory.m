function [x_data,y_data,t_data,eta_data,ui,vi] = b4_trajectory(p)

t = p.theta/(4*pi);
x_data = zeros(p.nimpacts,1); y_data = zeros(p.nimpacts,1); t_data = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 

xi=p.xi; yi=p.yi; ui=p.ui; vi=p.vi;

prefactor=  -p.d0*p.M(1)*p.G/(2*pi);

H_vec=zeros(p.Nx*p.Ny,1);
dH_vec= ones(p.Nx*p.Ny,1);

for n=1:p.nimpacts
    
    disp(['Impact number: ' num2str(n)])
   


    for nn=1:p.nsteps_impact 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec,dH_vec, t, p);

        

        t =t+ p.dt;
        if mod(nn,1)==0
            ploty = p.y; plotx = p.x(p.Nx/2);
            eta_plot  = zeros(size(ploty)); 
            for index=1:length(ploty)
                eta_plot(index) = sum(prefactor *p.K3_vec.* H_vec.*besselj(0,p.K_vec.* sqrt((plotx-xi).^2 +(ploty(index)-yi).^2 )));
            end
            plot(ploty,eta_plot);


            title(['Impact ' num2str(n) ' Step ' num2str(nn)]);
            ylim([-0.02,0.02])
            hold off;
            drawnow;
            pause(1/6);
        end
    end

end
end