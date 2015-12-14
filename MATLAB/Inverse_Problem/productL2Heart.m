function [ value ] = productL2Heart( u , v , MESH, FE_SPACE )
%PRODUCTL2HEART Computes the L2 scalar product in the heart
%   The reaction matrix must be assembled within the FE_SPACE data
%   structure! 

ubar = extend_with_zero( u , MESH ) ;
vbar = extend_with_zero( v , MESH ) ;

value = ubar' * FE_SPACE.A_reaction_heart * vbar ; 

end

