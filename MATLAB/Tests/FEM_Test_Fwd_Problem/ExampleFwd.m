%EXAMPLE2D shows how to solve a Pure Neumann problem with the FEM library.

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')

% Define which elements to use
fem = 'P2';

%% Import mesh

% meshFileName = '../../Mesh/Square/square_coarse' ;
% meshFileName = '../../Mesh/Square/square' ;
meshFileName = '../../Mesh/Square/square_fine' ;

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
data_file = 'DirichletNeumann_data' ;

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

% Print some information
fprintf(' * Number of Nodes           = %d \n',MESH.numNodes);
fprintf(' * Number of inner Nodes     = %d \n',MESH.numInnerNodes);

%% Define a control function 

ExampleLS ; % import circular level-set function

tau = 0.2 ;
w = LS( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) ) ;
w = 1 - smoothLS(w , tau) ;

% Visualize w
H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w' ) ; 
[X,Y] = meshgrid(-1:0.02:1) ; 
figure
surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp ; colormap jet ; title('Control function') ; axis equal ;



%% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

%  Visualize wbar
% figure
% pdeplot(MESH.vertices,[],MESH.elements(1:3,:) ,'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
%         'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
%         'colorbar','on','mesh','off');
% title('Extended control function')
% colormap(jet);

% USEFUL REMARK : the first MESH.numVertices elements of MESH.nodes
% correspond to MESH.vertices !!! 


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

%% Visualize diffusion coefficient

% [X,Y] = meshgrid( -2:0.1:2 , -2:0.1:2 ) ; 
% figure
% surf( X , Y , DATA.diffusion(X,Y) ) 
% shading interp
% colormap jet
% colorbar
% title('Diffusion coefficient')
% clear X; clear Y;


%% Assemble stiffness matrix
fprintf('\n Assembling ... ');
t_assembly = tic;
A              =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

%% Assemble rhs matrix 
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
A_source             =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);

%% Evaluate rhs 

F_source = DATA.coeffRhs * A_source * wbar ;

% Visualize source term 
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',F_source(1:MESH.numVertices),'xystyle','interp',...
       'zdata',F_source(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off' );
colormap(jet);
lighting phong

%% Apply boundary conditions

fprintf('\n Apply boundary conditions ');

[A_in, F_in, u_D]   =  ApplyBC_2D(A, F_source, FE_SPACE, MESH, DATA);

%% Impose zero-mean condition
A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
B = A_react * ones( MESH.numNodes , 1 ) ; 

A_total = [ A_in , B ; B' , 0 ] ; 
F_total = [ F_in ; 0 ] ; 

%% Solve
fprintf('\n Solve Au = f ... ');
t_solve = tic;
u                         = zeros(MESH.numNodes,1);
u_total = A_total \ F_total ; 
u(MESH.internal_dof)      = u_total(1 : end -1 ) ;
u(MESH.Dirichlet_dof)     = u_D;t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

%% Visualize solution 

figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',u(1:MESH.numVertices),'xystyle','interp',...
       'zdata',u(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
colormap(jet);
lighting phong
% axis equal




