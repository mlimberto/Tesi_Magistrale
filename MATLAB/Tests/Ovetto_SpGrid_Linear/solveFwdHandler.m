function u = solveFwdHandler( MESH , FE_SPACE , DATA , w , param , A_lhs , B , A_rhs )
%SOLVEFWDHANDLER Handler for the solveFwd function 
%   This function is responsible for the modification of the DATA 
%   of the problem according to the parameters

% Initialize stuff if needed
if nargin < 6 
    A_lhs = [] ; 
end
if nargin < 7
    B = [] ;
end
if nargin < 8 
    A_rhs = [] ;
end

% Modify DATA according to parameter

DATA.M0 = param(1) ; 
DATA.Mi = param(2) ;

DATA.diffusion = @(x,y,t,param)( DATA.M0 + (DATA.Mi + DATA.Me - DATA.M0)*( 1 - smoothLS( DATA.heartLS(x,y) , DATA.tauDiff) ) + 0.*x.*y);

DATA.coeffRhs = -1. * (DATA.vTr_i - DATA.vTr_e )*DATA.Mi ;

% Call the solveFwd function
u = solveFwd( MESH , FE_SPACE , DATA , w , A_lhs , B , A_rhs ) ;

end