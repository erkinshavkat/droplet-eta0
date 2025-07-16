function hole_IF(xi,yi,ui,vi,gamma,H,eta0,Nx,Lx)
%% Drop's initial position
disp('--- setting up ---')
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,0);
p.xi = xi; p.yi = yi; p.ui= ui; p.vi = vi;
p.nimpacts = 100;      % Number of impacts %%%%% ONLY THING TO EDIT HERE
lambdaf = p.lambdaf;
mem = p.mem;
%% Execute
fname = [pwd,'/results/output_double',num2str(Nx+Lx,'%.1f'),'mm.mat'];
[x_data,y_data,t_data, eta_data]  = compare_wave(p); 
eta_data = mean(eta_data,3);
% Save data
save(fname,'eta_data','x_data','y_data','t_data','gamma','H','Nx','Lx','mem');
end




