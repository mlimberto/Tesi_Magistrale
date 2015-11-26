%EXAMPLE2D shows how to solve a Pure Neumann problem with the FEM library.

% clc
% clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Define which elements to use
fem = 'P1';

%% Import mesh

meshFileName = '../../Mesh/Square/square_coarse' ;
% meshFileName = '../../Mesh/Square/square' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

% Add flags to identify mesh elements (in accordance with .geo file)

FLAG_HEART_REGION = 10 ;
FLAG_TORSO_REGION = 11 ;

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


%% Assemble matrix and right-hand side
fprintf('\n Assembling ... ');
t_assembly = tic;
[A, F]              =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);


%% Apply boundary conditions
fprintf('\n Apply boundary conditions ');

[A_in, F_in, u_D]   =  ApplyBC_2D(A, F, FE_SPACE, MESH, DATA);

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
       'colorbar','on','mesh','on');
colormap(jet);
lighting phong
axis equal

%% Compute L2 and H1 error

[errorL2,errorH1] = error2D(u, MESH, DATA, FE_SPACE);
    fprintf(' L2-error : %e H1-error : %e\n',errorL2, errorH1);




