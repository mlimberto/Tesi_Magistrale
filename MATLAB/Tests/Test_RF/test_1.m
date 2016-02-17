%% Import libraries

addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 
addpath( genpath( '../../RF' ) ) 

%% General parameters

NS = 10 ; % number of samples

N = 2 ; % number of KL bases

%% Define Gaussian covariance structure

corr.name = 'gauss';
corr.c0 = 1 ; 
corr.sigma = 0.1; 

%% Physical domain

mesh = linspace(-1,1,101)' ;   % generate a mesh
data.x = [-1; 1]; data.fx = [0; -1];    % specify boundaries

%% Generate random field information and samples

[F,KL] = randomfield( corr , mesh , 'nsamples' , NS , 'trunc' , N ,  'data' , data , 'filter' , 0.95 ) ;

% Plot bases

figure(2)
for i=1:N
   plot( mesh , KL.bases(:,i) ) 
   hold on 
end


%% Plot samples 
figure(1)
for i=1:NS
   plot( mesh , F(:,i) );
   hold on 
end
title('Field samples')

%% Generate more samples
 
W = randn(N,NS); 
F2 = repmat(KL.mean,1,NS) + KL.bases*diag(KL.sv)*W;

figure(1)
for i=1:NS
   plot( mesh , F(:,i) );
   hold on 
end

%% Create a sparse grid 

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

figure
for i=1:NS
   plot( mesh , F(:,i) );
   hold on 
end

