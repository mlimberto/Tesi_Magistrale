function [ u ] = solveFwd( MESH , FE_SPACE , DATA , w , A_lhs , B , A_rhs )
%SOLVEFWD Solves the simplified forward problem using a control function w
%   solveFwd( MESH , FE_SPACE , DATA , w , A_lhs , B , A_rhs )
%   The last three arguments are optional : the matrices and vectors
%   will be assembled if they are not provided by the user but this 
%   can take some additional time.

% Assemble stuff if needed
if nargin< 5 || isempty(A_lhs)
    A_lhs         =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
end
if nargin< 6 || isempty(B)
    A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
    B = A_react * ones( MESH.numNodes , 1 ) ; 
    clear A_react ;
end
if nargin< 7 || isempty(A_rhs)
    % Temporary replace the diffusion coefficient 
    temp = DATA.diffusion ;
    DATA.diffusion = @(x,y,t,param) 1 + 0*x.*y ;
    % Assemble matrix
    A_rhs         =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
    % Restore original diffusion coefficient
    DATA.diffusion = temp ;
end

% Extend control function to outer boundary
wbar = extend_with_zero( w , MESH ) ; 

% Evaluate rhs 
F = DATA.coeffRhs * A_rhs * wbar ;

% Apply boundary conditions
% There is no need to apply boundary conditions for this problem!!
% [A_in, F_in, u_D]   =  ApplyBC_2D(A_lhs, F, FE_SPACE, MESH, DATA);
% Instead we use A_lhs rather than A_in
%                F     rather than F_in 


% Impose zero-mean condition
A_total = [ A_lhs , B ; B' , 0 ] ; 
F_total = [ F ; 0 ] ; 

% Solve
u                         = zeros(MESH.numNodes,1);
u_total = A_total \ F_total ; 
u(MESH.internal_dof)      = u_total(1 : end -1 ) ;
% u(MESH.Dirichlet_dof) = u_D ; 
clear u_total;

end

