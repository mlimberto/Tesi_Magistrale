function [ wbar ] = extend_with_zero( w , MESH )
%EXTEND_WITH_ZERO provide a zero-extension on the outer mesh of a function defined in the inner mesh.

wbar = zeros( MESH.numNodes , 1 ) ;
wbar( MESH.indexInnerNodes , : ) = w ;  

end

