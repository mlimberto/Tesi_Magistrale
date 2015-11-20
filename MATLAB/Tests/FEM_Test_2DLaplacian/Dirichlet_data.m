%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique F??d??rale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Time 

data.dt = 0.01;
data.T = 1;
data.theta = 0.5;

% Source term
data.force = @(x,y,t,param)(8*pi^2*(sin(2*pi*x).*sin(2*pi*y)));

% Dirichlet
data.bcDir = @(x,y,t,param)(sin(2*pi*x).*sin(2*pi*y) + 0*x.*y);  

% Neumann
data.bcNeu = @(x,y,t,param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [1 2 3 4];
data.flag_neumann   = [];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x,y,t,param)(1 + 0.*x.*y);

% transport vector (first and second components)
data.transport{1} = @(x,y,t,param)(0 + 0.*x.*y);
data.transport{2} = @(x,y,t,param)(0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,t,param)(0 + 0.*x.*y);


% exact solution
data.uexact         = @(x,y,t,param)( sin(2*pi*x).*sin(2*pi*y));
data.uxexact        = @(x,y,t,param)( 2*pi*(cos(2*pi*x).*sin(2*pi*y)));
data.uyexact        = @(x,y,t,param)( 2*pi*(sin(2*pi*x).*cos(2*pi*y)));