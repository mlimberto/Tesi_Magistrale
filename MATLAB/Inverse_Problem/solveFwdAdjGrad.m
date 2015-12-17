function [ J , dw , u , p , normgradL2 , normgradH1 ] = solveFwdAdjGrad( w , MESH , FE_SPACE , DATA , zd )

% Assemble stuff if needed
if isempty(FE_SPACE.A_diffusion_total)
    disp 'Assembling total diffusion matrix'
    FE_SPACE.A_diffusion_total     =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
end
if isempty(FE_SPACE.B)
    disp 'Assembling B vector'
    A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
    FE_SPACE.B = A_react * ones( MESH.numNodes , 1 ) ; 
    clear A_react ;
end
if isempty(FE_SPACE.A_diffusion_heart)
    disp 'Assembling heart diffusion matrix'
    FE_SPACE.A_diffusion_heart    =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
end
if isempty(FE_SPACE.A_reaction_heart)
   disp 'Assembling heart reaction matrix' 
   FE_SPACE.A_reaction_heart = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' , [] , [] , DATA.FLAG_HEART_REGION ) ; 

end

% Solve forward problem 
wbar = extend_with_zero( w , MESH ) ; 
F = DATA.coeffRhs * FE_SPACE.A_diffusion_heart * wbar ;

A_total = [ FE_SPACE.A_diffusion_total , FE_SPACE.B ; FE_SPACE.B' , 0 ] ; 
F_total = [ F ; 0 ] ; 

u                         = zeros(MESH.numNodes,1);
u_total = A_total \ F_total ; 
u(MESH.internal_dof)      = u_total(1 : end -1 ) ;
clear u_total;


% Solve adjoint problem
F_adj = Apply_AdjBC( FE_SPACE , MESH , zd - u ) ;
F_total = [ F_adj ; 0 ] ; 

p                         = zeros(MESH.numNodes,1);
p_total = A_total' \ F_total ; 
p(MESH.internal_dof)      = p_total(1 : end -1 ) ;

% Evaluate objective function

J = eval_ObjFunction(MESH , DATA , FE_SPACE , w , u , zd , -F_adj ) ;


% Evaluate gradient
A_grad = FE_SPACE.A_diffusion_heart + FE_SPACE.A_reaction_heart ;

if nargout > 1
    F_grad = - DATA.coeffRhs * FE_SPACE.A_diffusion_heart * p ...
           + DATA.betaL2 * FE_SPACE.A_reaction_heart * wbar ...
           + DATA.betaGr * FE_SPACE.A_diffusion_heart * wbar ;
       
    dw = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ F_grad(MESH.indexInnerNodes) ;
   
    normgradL2 = sqrt( productL2Heart( dw , dw , MESH , FE_SPACE ) ) ;
    normgradH1 = sqrt( productH1Heart( dw , dw , MESH , FE_SPACE ) ) ;

end


end