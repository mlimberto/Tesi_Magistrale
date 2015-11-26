%EXAMPLESOURCETERM shows how to define the source term with \int{ \nabla w \nabla \phi }

% clc
% clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Define which elements to use
fem = 'P1';

%% Import mesh

% meshFileName = '../../Mesh/Square/square_coarse' ;
meshFileName = '../../Mesh/Square/square' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

% Add flags to identify mesh elements (in accordance with .geo file)

FLAG_HEART_REGION = 10 ;
FLAG_TORSO_REGION = 11 ;

FLAG_TORSO_BOUNDARY_DIRI = 1 ; 
FLAG_TORSO_BOUNDARY_NEU = 2 ;
FLAG_HEART_BOUNDARY = 3 ;


%% Visualize mesh

%  figure
%  pdeplot(vertices,[],elements(1:3,:))
%  axis equal

% figure
% pdeplot(vertices , [] , elements(1:3 ,elements(4,:)==FLAG_HEART_REGION ) ) ;
% title('Heart domain') ;
% 
% figure
% pdeplot(vertices , [] , elements(1:3 ,elements(4,:)==FLAG_TORSO_REGION ) ) ;
% title('Torso domain') ;


%% Set parameters from file
data_file = 'Dirichlet_data' ;

DATA   = read_DataFile(data_file);
DATA.param = [] ;
t = [] ;

%% Fill MESH data structure
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements  = elements;
MESH.numVertices = size(vertices,2);

% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(strcat('P1to',fem,'mesh','2D'),elements, vertices, boundaries);
else
    MESH.nodes = vertices;
end

%% Update Mesh data with BC information and geometrical maps
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

%% Fill inner mesh data structure

MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== FLAG_HEART_REGION ) ; 
MESH.indexInnerBoun = find( MESH.boundaries( 5 , :) == FLAG_HEART_BOUNDARY ) ;
MESH.indexInnerVert = MESH.boundaries(1, MESH.indexInnerBoun ) ;

MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;

MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 

MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% %Verify that the internal nodes have been properly identified 
% figure
% plot( MESH.nodes(1,:) , MESH.nodes(2,:) , 'o' )
% hold on 
% plot( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) , '*' ) 
% legend('Nodes' , 'Inner Nodes')

%% Create and fill the FE_SPACE data structure
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

%% Define a control function and extend it to the whole mesh

% Use a circular level set function
R = 0.25 ;
LS = @(x,y)( sqrt( x.*x + y.*y) - R ) ;

tau = 0.4 ;
w = LS( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) ) ;
w = 1 - smoothLS(w , tau) ;

wbar = extend_with_zero( w , MESH) ; 

%  Visualize wbar
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:) ,'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
        'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
        'colorbar','on','mesh','on');
title('Extended control function')
colormap(jet);


%% Assemble stiffness matrix 
fprintf('\n Assembling stiffness matrix ... ');
t_assembly = tic;
A            =  Assembler_2D(MESH, DATA, FE_SPACE);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

%% Assemble rhs matrix 
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
A_source             =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [1 1] , [] );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);

% REMARK Since w_bar is zero on the torso, we could just assemble the
% matrix on the heart subdomain using the expression
% Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [1 1] , [] ,
% FLAG_HEART_REGION ) ;

% Check that the source term in this case is the same
A_test = Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [1 1] , [] , FLAG_HEART_REGION ) ;
sourceError=  norm( A_test * wbar - A_source * wbar ) ;
assert( sourceError < 1e-14 );
clear A_test

%% Evaluate rhs 

F_source = A_source * wbar ;

% Visualize source term
figure
Ff = full(F_source) ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',Ff(1:MESH.numVertices),'xystyle','interp',...
       'zdata',Ff(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on','mesh','on');
colormap(jet);
lighting phong
title('Source term')
clear Ff

%% Apply boundary conditions
fprintf('\n Apply boundary conditions ');

[A_in, F_in, u_D]   =  ApplyBC_2D(A, F_source, FE_SPACE, MESH, DATA);

%% Solve
fprintf('\n Solve Au = f ... ');
t_solve = tic;
u                         = zeros(MESH.numNodes,1);
u(MESH.internal_dof)      = A_in \ F_in;
u(MESH.Dirichlet_dof)     = u_D;t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

%% Visualize solution 

figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',u(1:MESH.numVertices),'xystyle','interp',...
       'zdata',u(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on','mesh','on');
colormap(jet);
title('Solution')
lighting phong



