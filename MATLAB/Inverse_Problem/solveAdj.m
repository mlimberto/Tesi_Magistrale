function [ p , F_adj ] = solveAdj( MESH , FE_SPACE , DATA , w , u , zd, A_lhs , B  )
%SOLVAADJ solves the simplified adjoint problem using a control function w 
%   solveAdj( MESH , FE_SPACE , DATA , w , u , zd , A_lhs , B , A_rhs )
%   The last three arguments are optional : the matrices and vectors
%   will be assembled if they are not provided by the user but this 
%   can take some additional time.


% Assemble stuff if needed
if nargin< 7 || isempty(A_lhs)
    A_lhs         =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
    A_lhs         =  A_lhs' ;
end
if nargin< 8 || isempty(B)
    A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
    B = A_react * ones( MESH.numNodes , 1 ) ; 
    clear A_react ;
end

% Extend control function to outer boundary
% wbar = extend_with_zero( w , MESH ) ; 

% Get source term
F_adj = Apply_AdjBC( FE_SPACE , MESH , zd - u ) ;

% Impose zero-mean condition
A_total = [ A_lhs , B ; B' , 0 ] ; 
F_total = [ F_adj ; 0 ] ; 

% Solve
p                         = zeros(MESH.numNodes,1);
p_total = A_total \ F_total ; 
p(MESH.internal_dof)      = p_total(1 : end -1 ) ;

end

