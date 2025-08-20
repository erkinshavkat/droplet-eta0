function hole_IF(xi,yi,ui,vi,gamma,H,Nx,Lx,Nk,kmin,kmax)
%% Drop's initial position
disp('--- setting up ---')
eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,kmin,kmax);
p.xi = xi; p.yi = yi; p.ui= ui; p.vi = vi;
p.nimpacts = 30;      % Number of impacts %%%%% ONLY THING TO EDIT HERE
lambdaf = p.lambdaf;
mem = p.mem;
%% Execute

resultdir= 'results/b4vFaria_1bc/';
mkdir(resultdir)
fname = [resultdir,sprintf('faria_N%d_Nk%dx_Lk%dx.mat',Nx,kmax)];
[eta_b4, eta_faria]  = compare_wave(p); 
%save(fname,'eta_b4','eta_faria','gamma','H','Nx','Lx','mem');
end




