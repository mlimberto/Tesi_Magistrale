%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,t,param)(0.5*pi^2*(sin(0.5*pi*x).*sin(0.5*pi*y)));
% data.force = @(x,y,t,param)( 0*x.*y) ;

% Dirichlet
data.bcDir = @(x,y,t,param)(sin(0.5*pi*x).*sin(0.5*pi*y) + 0*x.*y);  
% data.bcDir = @(x,y,t,param)( x +  y );

% Neumann
data.bcNeu = @(x,y,t,param)(- 2*pi*(cos(2*pi*x).*sin(2*pi*y)) + 0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)( 0.*x.*y);
data.bcRob_fun      = @(x,y,t,param)( 0.*x.*y);

% BC flag
data.flag_dirichlet = [ 1 2 ]; % This is a pure Dirichlet problem
data.flag_neumann   = [ ];
data.flag_robin     = [ ];

% diffusion
data.diffusion = @(x,y,t,param)(1 + 0.*x.*y);

% transport vector (first and second components)
data.trasport{1} = @(x,y,t,param)(0 + 0.*x.*y);
data.trasport{2} = @(x,y,t,param)(0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,t,param)(0 + 0.*x.*y);


% exact solution    = @(x,y,t,param)( 1 + 0*x.*y );
data.uexact         = @(x,y,t,param)( sin(0.5*pi*x).*sin(0.5*pi*y));
data.uxexact        = @(x,y,t,param)( 0.5*pi*(cos(0.5*pi*x).*sin(0.5*pi*y)));
data.uyexact        = @(x,y,t,param)( 0.5*pi*(sin(0.5*pi*x).*cos(0.5*pi*y)));


