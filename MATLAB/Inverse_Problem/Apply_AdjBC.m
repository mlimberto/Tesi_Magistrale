function [ F_in ] = Apply_AdjBC( FE_SPACE , MESH , f  )
%APPLY_ADJBC Applies the Neumann boundary condition to the adjoint problem
%matrix
%   We use it to evaluate the rhs term (boundary integral of u - zd)
%   f has to be replaced with -1 * (u-zd )

% Evaluate quadrature nodes, weights and basis functions
[csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
[phi]          =  fem_basis2D(FE_SPACE.fem, csi, 0*csi, 1);
eta            =  1 - csi;
nqn            =  length(csi);
    
nof         = length(MESH.Neumann_side);
nbn         = MESH.numBoundaryDof;
    
Rrows       = zeros(nbn*nof,1);
Rcoef       = Rrows;
    
xlt = zeros(nof,nqn); ylt = xlt;
coord_ref = [eta; csi];

for j = 1 : 2
    dof = MESH.boundaries(j,MESH.Neumann_side);
    vtemp = MESH.vertices(1,dof);
    xlt = xlt + vtemp'*coord_ref(j,:);
    vtemp = MESH.vertices(2,dof);
    ylt = ylt + vtemp'*coord_ref(j,:);
end

dataNeu = scatteredInterpolant( MESH.nodes(1,:)' , MESH.nodes(2,:)' , f ) ; % verify the sign!

u_Neumann = dataNeu(xlt,ylt);
one       = ones(nof,nqn);
u_Neumann = u_Neumann.*one;

x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Neumann_side));
y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Neumann_side));
   
side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
    
for l = 1 : nof
    face = MESH.Neumann_side(l);
        
    u_Neumann_loc  = u_Neumann(l,:).*wi;
    u_Neumann_loc  = u_Neumann_loc(1,:)';
      
    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
    Rcoef(1+(l-1)*nbn:l*nbn)    = side_length(l)*phi*u_Neumann_loc;
end

F_in = sparse(Rrows,1,Rcoef,MESH.numNodes,1);

end

