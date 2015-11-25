function [ y ] = smoothLS( x , tau )
%SMOOTHLS is an approximation of a Heaviside function with parameter tau

y = ones(  size(x) ).*( x > tau ) + ...
    zeros( size(x) ).*( x < tau ) + ...
    0.5 .* ( 1 + x./tau + sin( pi.*x./tau)./pi ).*( abs(x) <= tau) ; 

end

