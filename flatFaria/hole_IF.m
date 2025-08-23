function hole_IF(xi,yi,ui,vi,gamma,H,Nx,Lx,Nk,kmin,kmax,theta)
%% Drop's initial position
disp('--- setting up ---')
eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,kmin,kmax,theta);
p.xi = xi; p.yi = yi; p.ui= ui; p.vi = vi;
p.nimpacts = 20;      % Number of impacts %%%%% ONLY THING TO EDIT HERE
lambdaf = p.lambdaf;
mem = p.mem;
%% Execute

resultdir= 'results/varparam/b1/';
mkdir(resultdir)
fname = [resultdir,sprintf('b1_N%d_L%d_H%d_mem%.2f_ph%.2f.mat',Nx,Lx,1000*H,mem,theta)];
[eta_data, eta_intermediate]  = compare_H(p);
%save(fname,'x_data','y_data','eta_data','eta_intermediate','theta','gamma','H','Nx','Lx','mem');

end




