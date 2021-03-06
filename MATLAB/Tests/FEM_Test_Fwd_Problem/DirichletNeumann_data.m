%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Physiological and geometrical parameters

data.heartLS = @(x,y)( max( abs(x) , abs(y) ) - 1 ) ;

data.M0 = 2.39 ;
data.Mi = 3.0  ; 
data.Me = 2.0  ;

data.tildeMi = 3.0 ;
data.tildeMe = 1.0 ;

data.vTr_i = -60 ;  
data.vTr_e = -90 ;

data.coeffRhs = -1. * (data.vTr_i - data.vTr_e )*data.Mi ;

% Source term
% data.force = @(x,y,t,param)(0.5*pi^2*(sin(0.5*pi*x).*sin(0.5*pi*y)));
data.force = @(x,y,t,param)( 0*x.*y) ;

% Dirichlet
% data.bcDir = @(x,y,t,param)( x +  y );

% Neumann
% data.bcNeu = @(x,y,t,param)(1+ 0*x.*y );
data.bcNeu = @(x,y,t,param)( 0*x.*y );

% Robin
% data.bcRob_alpha    = @(x,y,t,param)( 0.*x.*y);
% data.bcRob_fun      = @(x,y,t,param)( 0.*x.*y);

% BC flag
data.flag_dirichlet = [];
data.flag_neumann   = [1 2];
data.flag_robin     = [ ];

% diffusion
data.tauDiff = 0.08 ;
data.diffusion = @(x,y,t,param)( data.M0 + (data.Mi + data.Me - data.M0)*( 1 - smoothLS( data.heartLS(x,y) , data.tauDiff) ) + 0.*x.*y);

% transport vector (first and second components)
data.trasport{1} = @(x,y,t,param)(0 + 0.*x.*y);
data.trasport{2} = @(x,y,t,param)(0 + 0.*x.*y);

% reaction
% data.reaction = @(x,y,t,param)(0 + 0.*x.*y);
data.reaction = @(x,y,t,param)(1 + 0.*x.*y); % this is used to impose the zero mean condition!!


% exact solution    = @(x,y,t,param)( 1 + 0*x.*y );
data.uexact         = @(x,y,t,param)( 0*x.*y );
data.uxexact        = @(x,y,t,param)( 0*x.*y );
data.uyexact        = @(x,y,t,param)( 0*x.*y );


