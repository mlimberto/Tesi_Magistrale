function [up] = solveFwdAdjHandler( MESH , FE_SPACE , DATA , w , zd , param , A_lhs_fwd , B , A_rhs_fwd )
%SOLVEFWDHANDLER Handler for the solveFwd and solveAdj functions. 
%   This function is responsible for the modification of the DATA 
%   of the problem according to the parameters
%   Then both the forward and adjoint problem are solved with 
%   respect to the chosen combination of parameters.

% Initialize stuff if needed
if nargin < 6 
    A_lhs_fwd = [] ; 
end
if nargin < 7
    B = [] ;
end
if nargin < 8 
    A_rhs_fwd = [] ;
end

% Modify DATA according to parameter

DATA.M0 = param(1) ; 
DATA.Mi = param(2) ;

DATA.coeffRhs = -1. * (DATA.vTr_i - DATA.vTr_e )*DATA.Mi ;

% Call the solveFwd function
u = solveFwd( MESH , FE_SPACE , DATA , w , A_lhs_fwd , B , A_rhs_fwd ) ;

% Call the solveAdj function
p = solveAdj( MESH , FE_SPACE , DATA , w , u , zd , [] , [] ) ;

% Merge the solutions together
up = [ u ; p ] ;

% REMARK : This could go much faster by pre-assembling the lhs and rhs matrices

end

