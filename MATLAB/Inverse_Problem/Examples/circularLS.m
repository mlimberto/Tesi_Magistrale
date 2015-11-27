function [ z ] = circularLS( x , y , R , D )
%CENTEREDCIRCULARLS Provides simple circular level set function

if (nargin < 3 || isempty(R) ) 
    R = 0.25 ; % set default radius
end

if (nargin < 4 || isempty(D) )
   D = [ 0 0 ] ; 
end

z = sqrt( (x-D(1) ).*(x-D(1) ) + (y-D(2)).*(y-D(2)) ) - R ;

end

