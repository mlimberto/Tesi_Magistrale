function out = solveFwdAdjHandler( MESH , FE_SPACE , DATA , w , zd , param  )
%SOLVEFWDHANDLER Handler for the solveFwd and solveAdj functions. 
%   This function is responsible for the modification of the DATA 
%   of the problem according to the parameters
%   Then both the forward and adjoint problem are solved with 
%   respect to the chosen combination of parameters.

% Modify DATA according to parameter
DATA.M0 = param(1) ; 
DATA.Mi = param(2) ;
DATA.diffusion = @(x,y,t,param)( DATA.M0 + (DATA.Mi + DATA.Me - DATA.M0)*( 1 - smoothLS( DATA.heartLS(x,y) , DATA.tauDiff) ) + 0.*x.*y);
DATA.coeffRhs = -1. * (DATA.vTr_i - DATA.vTr_e )*DATA.Mi ;

% Assemble stuff if needed

if isempty(FE_SPACE.B)
    B = [] ;
else
    B = FE_SPACE.B ;
end

if isempty(FE_SPACE.A_diffusion_heart)
    F_rhs_fwd = [] ;
    FE_SPACE.A_diffusion_heart = Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
else
    wbar = extend_with_zero( w , MESH ) ; 
    F_rhs_fwd = DATA.coeffRhs * FE_SPACE.A_diffusion_heart * wbar ;
end

% Assemble fwd lhs matrix
A_lhs_fwd         =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');



% Call the solveFwd function
u = solveFwd( MESH , FE_SPACE , DATA , w , A_lhs_fwd , B , F_rhs_fwd ) ;

% Call the solveAdj function
p = solveAdj( MESH , FE_SPACE , DATA , w , u , zd , A_lhs_fwd' , B ) ;

% Compute the gradient contribution term
F_grad = DATA.Mi .* (DATA.vTr_i - DATA.vTr_e ) * FE_SPACE.A_diffusion_heart * p ;

% Merge the solutions together
out = [F_grad ; u ; p ] ;

end
