%% Import libraries

addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 
addpath( genpath( '../../RF' ) ) 

%% General parameters

NS = 10 ; % number of samples

N = 2 ; % number of KL bases

%% Physical domain

N_grid = 11;

x = linspace(-1,1,N_grid);
[X,Y] = meshgrid(x,x); mesh = [X(:) Y(:)]; % 2-D mesh

%% Define Gaussian covariance structure

corr.name = 'gauss';
corr.c0 = [ 0.2 1] ; 
corr.sigma = cos(pi*mesh(:,1)).*sin(2*pi*mesh(:,2))+1.5;

%% Generate RF expansion

[F,KL] = randomfield(corr,mesh,'trunc', N);

% plot the realization
surf(X,Y,reshape(F,N_grid,N_grid)); colorbar;
shading interp

%% Generate gaussian sparse grid 

mu = 0 ; 
sigma = 1 ;
knots = @(n) knots_gaussian( n , mu , sigma ) ; % gaussian quadrature nodes
w = 2 ; % level of the grid

S = smolyak_grid( N ,w,knots,@lev2knots_doubling);
Sr=reduce_sparse_grid(S);

% Plot grid points (first two dimensions)
if ( N > 1 )
    figure
    scatter( Sr.knots(1,:) , Sr.knots(2,:) , 'filled' )
    grid on 
end

% Plot grid points (first three dimensions)
if ( N > 2 )
    figure
    scatter3( Sr.knots(1,:) , Sr.knots(2,:) , Sr.knots(3,:) , 'filled')
    grid on 
end

%% Evaluate random field at grid nodes

NS = size(Sr.knots,2) ;

F = repmat(KL.mean,1,NS) + KL.bases*diag(KL.sv)*Sr.knots ;

%% Plot some of the random fields

V = [1 2 6 9] ;

figure
for i=1:4
   subplot(2,2,i)
   surf(X,Y,reshape(F(:,V(i)),N_grid,N_grid)); colorbar;
   shading interp    
end





