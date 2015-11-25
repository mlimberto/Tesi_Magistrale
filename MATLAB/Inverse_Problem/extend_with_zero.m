function [ wbar ] = extend_with_zero( w , MESH )
%extend_with_zero Extends a function defined in the inner mesh to the outer
%mesh with zero

wbar = zeros( MESH.numNodes , 1 ) ;
wbar( MESH.indexInnerNodes , : ) = w ;  

end

