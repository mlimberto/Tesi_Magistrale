function [ y ] = smoothDLS( x , tau )
%SMOOTHDLS Provides an approximation of the Dirac Delta 

y = zeros( size(x) ).*( abs(x) > tau ) + ...
    0.5 ./ tau .* ( 1 + cos( pi.*x./tau) ).*( abs(x) <= tau) ; 

end

