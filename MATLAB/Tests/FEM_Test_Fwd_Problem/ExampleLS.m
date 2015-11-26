% Set radius of ischemic region
R = 0.25 ;

% Define circular level-set
LS = @(x,y)( sqrt( x.*x + y.*y) - R ) ;

% Plot level-set and smoothed level-set

% x = linspace(-1,1,100);
% [X,Y] = meshgrid(x,x) ;
% Z = LS(X,Y) ;
% 
% tau = 0.2;
% smoothedZ = smoothLS(Z , tau);

% surf(X,Y,Z)
% % colormap jet
% surfc(X,Y,Z)
% shading interp
% title('Circular level set function')

% figure
% surf(X,Y,1 - smoothedZ )
% % colormap jet
% shading interp
% title('Ischemia location')

