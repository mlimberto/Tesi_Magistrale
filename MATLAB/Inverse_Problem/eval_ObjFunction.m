function [ J ] = eval_ObjFunction( MESH , DATA , FE_SPACE , w , u , zd , F_u , A_w )
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
if isempty( A_w ) 
   sprintf('Assembling A_w matrix ...') 
   A_react  = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' , [] , [] , DATA.FLAG_HEART_REGION ) ; 
   A_source =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
   A_w = A_react + A_source ; 
end

J = J + 0.5 * DATA.beta * wbar' * A_w * wbar ; 

end

