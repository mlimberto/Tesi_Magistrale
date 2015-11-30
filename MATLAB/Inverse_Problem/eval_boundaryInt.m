function [ I ] = eval_boundaryInt( MESH , FE_SPACE , f , selected_boundary  )
%EVAL_BOUNDARYINT evaluates the integral on the torso boundary of a certain
%function f
%   This function uses a scattered interpolant so there is a small
%   interpolation error when evaluating the function since linear 
%   interpolation is used


if nargin < 4
   selected_boundary = 1:size(MESH.boundaries , 2) ; 
end

% First we evaluate quadrature nodes at the boundary
[csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
eta            =  1 - csi;
nqn            =  length(csi);
    
nof         = length(selected_boundary);

I = 0.0 ;
    
xlt = zeros(nof,nqn); ylt = xlt;
coord_ref = [eta; csi];
for j = 1 : 2
    dof = MESH.boundaries(j,selected_boundary);
    vtemp = MESH.vertices(1,dof);
    xlt = xlt + vtemp'*coord_ref(j,:);
    vtemp = MESH.vertices(2,dof);
    ylt = ylt + vtemp'*coord_ref(j,:);
end
    
% Then we evaluate the function we wish to integrate at quadrature nodes
H = scatteredInterpolant( MESH.nodes(1,:)' , MESH.nodes(2,:)' , f) ;
f_value = H(xlt,ylt);
one       = ones(nof,nqn);
f_value = f_value.*one;
    
x    =  MESH.vertices(1,MESH.boundaries(1:2, selected_boundary));
y    =  MESH.vertices(2,MESH.boundaries(1:2, selected_boundary));
    
    side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
    
    for l = 1 : nof
        
        u_Neumann_loc  = f_value(l,:).*wi;
        u_Neumann_loc  = u_Neumann_loc(1,:)';
        I    = I + side_length(l)*ones( 1, FE_SPACE.quad_order )*u_Neumann_loc ;
    end

end

