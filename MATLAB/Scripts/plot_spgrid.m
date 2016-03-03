% Plot full vs sparse grid 

clc
clear all
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 

%% 

N=2; % approximation of two variables
knots=@(n) knots_CC(n,-1,1,'nonprob'); % knots
w = 2; %level
S = smolyak_grid(N,w,knots,@lev2knots_doubling); % grid


[X,Y] = meshgrid(knots(5),knots(5));
X = X(:) ; Y = Y(:) ;

%%

h = figure
s1 = subplot(121)
plot(X,Y,'color','b','marker','o','MarkerFaceColor','b','Linestyle','None');
grid on 
s2 = subplot(122)
plot_grid(S,[],'color','b','marker','o','MarkerFaceColor','b');

%% 

s1.Position = [0.05 0.05 0.4 0.9] ;
s2.Position = [0.55 0.05 0.4 0.9] ;