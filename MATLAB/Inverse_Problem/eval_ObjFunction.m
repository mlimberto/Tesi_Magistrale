function [ J ] = eval_ObjFunction( MESH , DATA , FE_SPACE , w , u , zd , F_u )
%EVAL_OBJFUNCTION Evaluates the objective function J(u,w)

wbar = extend_with_zero( w , MESH ) ;

J = 0.0 ;

% Evaluate term dependent on u 

if isempty(F_u) 
   % We use the function that applies the Neumann boundary conditions
   % to evaluate the boundary integral
   F_u = Apply_AdjBC ( FE_SPACE , MESH , u - zd ) ; 
end

J = J + 0.5 * ( u - zd )' * F_u ;

% Evaluate term dependent on w 

if isempty( FE_SPACE.A_reaction_heart ) 
       FE_SPACE.A_reaction_heart  = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' , [] , [] , DATA.FLAG_HEART_REGION ) ; 
end

if isempty( FE_SPACE.A_diffusion_heart ) 
       FE_SPACE.A_diffusion_heart =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
end

J = J + 0.5 * DATA.betaL2 * wbar' * FE_SPACE.A_reaction_heart * wbar ...
      + 0.5 * DATA.betaGr * wbar' * FE_SPACE.A_diffusion_heart * wbar ; 

end

