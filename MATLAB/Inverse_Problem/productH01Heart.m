function [ value ] = productH01Heart( u , v , MESH, FE_SPACE )
%PRODUCTH01HEART Computes the H01 scalar product in the heart
%   The diffusion matrix must be assembled within the FE_SPACE data
%   structure! 

ubar = extend_with_zero( u , MESH ) ;
vbar = extend_with_zero( v , MESH ) ;

value = ubar' * FE_SPACE.A_diffusion_heart * vbar ; 

end

