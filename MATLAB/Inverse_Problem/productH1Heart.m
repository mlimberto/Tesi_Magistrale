function [ value ] = productH1Heart( u , v , MESH, FE_SPACE )
%PRODUCTH1HEART Computes the H01 scalar product in the heart
%   The diffusion and reaction matrices must be assembled within the FE_SPACE data
%   structure! 

value =  productL2Heart( u , v , MESH , FE_SPACE ) ...
      + productH01Heart( u , v , MESH , FE_SPACE ) ; 

end

