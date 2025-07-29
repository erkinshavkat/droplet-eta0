function hole_IF(xi,yi,ui,vi,gamma,H,eta0,Nx,Lx,Nk_multiple)
%% Drop's initial position
disp('--- setting up ---')
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk_multiple);
p.xi = xi; p.yi = yi; p.ui= ui; p.vi = vi;
p.nimpacts = 10;      % Number of impacts %%%%% ONLY THING TO EDIT HERE
lambdaf = p.lambdaf;
mem = p.mem;
%% Execute

resultdir= 'results/b4vFaria_1bc/';
mkdir(resultdir)
fname = [resultdir,sprintf('faria_N%d_Nk%dx.mat',Nx,Nk_multiple)];
[x_data,y_data,t_data, eta_data]  = trajectory_IF_matt(p); 
save(fname,'eta_data','x_data','y_data','t_data','gamma','H','Nx','Lx','mem');

fname = [resultdir,sprintf('b4_N%d_Nk%dx.mat',Nx,Nk_multiple)];
[x_data,y_data,t_data, eta_data]  = b4_trajectory(p); 
save(fname,'eta_data','x_data','y_data','t_data','gamma','H','Nx','Lx','mem');
end




