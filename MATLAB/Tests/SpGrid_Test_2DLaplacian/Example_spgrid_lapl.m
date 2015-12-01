%EXAMPLE_SPGRID_LAPL shows how to solve a laplacian with the FEM library 
% but the location of the source term is random (uniformly sampled)

% clc
% clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Add path to Sparse grid Matlab toolkit
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 

% Define which elements to use
fem = 'P1';

%% Physical domain

% meshFileName = '../../Mesh/Square/square_coarse' ;
meshFileName = '../../Mesh/Square/square' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

% pdeplot(vertices,[],elements(1:3,:))
% axis equal

DATA   = read_DataFile('DirichletNeumann_data');
DATA.param = [] ; 
DATA.t = [];

MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements  = elements;
MESH.numVertices = size(vertices,2);

if ~strcmp(fem,'P1')
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(strcat('P1to',fem,'mesh','2D'),elements, vertices, boundaries);
else
    MESH.nodes = vertices;
end

% Update Mesh data with BC information and geometrical maps
[numElemDof,numBoundaryDof]  = select2D(fem);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;

% Update MESH with BC information
[MESH]         = BC_info(MESH, DATA);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf2D(MESH.vertices, MESH.elements);   

% Compute quadrature nodes and weights on the reference element
quad_order                  = 4; 
[quad_nodes, quad_weights]  = dunavant_quad(quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis2D('P1', quad_nodes(1,:), quad_nodes(2,:));

% Create and fill the FE_SPACE data structure
FE_SPACE.fem              = fem;
FE_SPACE.numDof           = length(MESH.internal_dof);
FE_SPACE.numElemDof       = numElemDof;
FE_SPACE.numBoundaryDof   = numBoundaryDof;

% Store quadrature nodes and weights on the reference element
FE_SPACE.quad_order    = quad_order;
FE_SPACE.quad_nodes    = quad_nodes;
FE_SPACE.quad_weights  = quad_weights;
FE_SPACE.numQuadNodes  = length(FE_SPACE.quad_nodes);

% Evaluate basis functions in the quad points on the reference element
[FE_SPACE.phi, FE_SPACE.dcsiphi, FE_SPACE.detaphi]  =  ...
    fem_basis2D(FE_SPACE.fem, FE_SPACE.quad_nodes(1,:), FE_SPACE.quad_nodes(2,:));


fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf('-------------------------------------------\n');


%% Stochastic domain 

% The location of the source term is uniformly distributed in the 
%[-1 1]x[-1 1] domain. 

% We use a sparse grid to discretize this space

N=2; % approximation of two variables
knots=@(n) knots_CC(n,-1,1,'nonprob'); % knots
w = 2; %level

[ lev2knots , idxset ] = define_functions_for_rule( 'SM' , N ) ;

S = smolyak_grid(N,w,knots,lev2knots , idxset); % grid

% plot the grid itself
figure
plot_grid(S,[],'color','k','marker','o','MarkerFaceColor','k');

% Reduce grid
Sr = reduce_sparse_grid( S ) ;

%% Solve diffusion problem

Q = size(Sr.knots , 2 ) ;

U = [] ;

for i = 1:Q
    source_location = Sr.knots( : , i) ;
    
    u = solveDiffusion( MESH  , FE_SPACE , DATA , source_location ) ;
    
    U = [ U u ] ;
    
end

figure
for i  =1:3
    for j = 1:3
        subplot(3,3, (i-1)*3 + j )
        u = U(: , 3*(i-1)+j ) ;
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',u(1:MESH.numVertices),'xystyle','interp',...
            'zdata',u(1:MESH.numVertices),'zstyle','continuous',...
            'colorbar','on','mesh','on');
        colormap(jet);
        lighting phong
    end
end



