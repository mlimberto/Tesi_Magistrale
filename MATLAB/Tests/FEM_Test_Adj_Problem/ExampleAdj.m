%EXAMPLEADJ solves the adjoint problem with the FEM library.

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First part : Mesh and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define which elements to use
fem = 'P1';

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

%% Assemble stiffness matrix
fprintf('\n Assembling ... ');
t_assembly = tic;
A              =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

%% Assemble rhs matrix 
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
A_source             =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , FLAG_HEART_REGION );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);

%% Build vector for zero-mean condition
A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
B = A_react * ones( MESH.numNodes , 1 ) ; 
clear A_react ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second part : Solve the fwd problem to find zd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Evaluating target solution zd ... \n');

%% Define a target control function 

tau = 0.2 ; % smoothing level
R = 0.4 ; % radius
D = [1.1 0]; % displacement

w = circularLS( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) , R , D ) ;
w = 1 - smoothLS(w , tau) ;

% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
% H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w' ) ; 
% [X,Y] = meshgrid(-1:0.02:1) ; 
% figure
% surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% shading interp ; colormap jet ; title('Target control function') ; axis equal ;
% clear H ; clear X ; clear Y ; 

%% Evaluate rhs 

F_source = DATA.coeffRhs * A_source * wbar ;

% Visualize rhs term 
% figure
% pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',F_source(1:MESH.numVertices),'xystyle','interp',...
%        'zdata',F_source(1:MESH.numVertices),'zstyle','continuous',...
%        'colorbar','on', 'mesh' , 'off' );
% colormap(jet);
% lighting phong
% title('Source term')

%% Apply boundary conditions

fprintf('\n Apply boundary conditions ');

[A_in, F_in, u_D]   =  ApplyBC_2D(A, F_source, FE_SPACE, MESH, DATA);

%% Impose zero-mean condition

A_total = [ A_in , B ; B' , 0 ] ; 
F_total = [ F_in ; 0 ] ; 

%% Solve
fprintf('\n Solve Au = f ... ');
t_solve = tic;
zd                         = zeros(MESH.numNodes,1);
zd_total = A_total \ F_total ; 
zd(MESH.internal_dof)      = zd_total(1 : end -1 ) ;
zd(MESH.Dirichlet_dof)     = u_D;t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);
clear zd_total;

%% Visualize solution 

figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
       'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
colormap(jet); lighting phong; title('Target solution z_d')
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third part : Solve the fwd problem and find u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Solving the state problem to find u ... \n');

%% Define a control function 

% clear previous data
clear w ; clear wbar ;

tau = 0.2 ; % smoothing level

% Since we use a level set we define radius and center coordinates
R = 0.4 ; % radius
D = [0 0]; % displacement

w = circularLS( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) , R , D ) ;
w = 1 - smoothLS(w , tau) ;

% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w' ) ; 
[X,Y] = meshgrid(-1:0.02:1) ; 
figure
surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp ; colormap jet ; title('Control function') ; axis equal ;
clear H ; clear X ; clear Y ; 

%  Visualize wbar
% figure
% pdeplot(MESH.vertices,[],MESH.elements(1:3,:) ,'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
%         'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
%         'colorbar','on','mesh','off');
% title('Extended control function')
% colormap(jet);

% USEFUL REMARK : the first MESH.numVertices elements of MESH.nodes
% correspond to MESH.vertices !!! 

% Visualize diffusion coefficient

% [X,Y] = meshgrid( -2:0.1:2 , -2:0.1:2 ) ; 
% figure
% surf( X , Y , DATA.diffusion(X,Y) ) 
% shading interp
% colormap jet
% colorbar
% title('Diffusion coefficient')
% clear X; clear Y;

%% Assemble stiffness matrix

% We use the stiffness matrix we assembled for z_d
% In the case of the non-linear problem this can't be done and 
% we would have to reassemble the stiffness matrix every time we 
% update the control function w

%% Evaluate rhs 

% clear previous data
clear A_in ; clear A_total ; clear F_source ; clear F_in ; clear F_total ;  

F_source = DATA.coeffRhs * A_source * wbar ;

% Visualize rhs term 
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',F_source(1:MESH.numVertices),'xystyle','interp',...
       'zdata',F_source(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off' );
colormap(jet);
lighting phong
title('Source term')

%% Apply boundary conditions

fprintf('\n Apply boundary conditions ');

[A_in, F_in, u_D]   =  ApplyBC_2D(A, F_source, FE_SPACE, MESH, DATA);

%% Impose zero-mean condition

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
colormap(jet); lighting phong ; title('Solution of state problem u') 
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth part : Solve the adj problem and find p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assemble adjoint problem stiffness matrix 

A_adj = A' ; 

%% Evaluate volume source term 

F_adj = zeros(size(A_adj,1),1);

%% Apply Neumann boundary conditions

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

dataNeu = scatteredInterpolant( MESH.nodes(1,:)' , MESH.nodes(2,:)' , zd -u ) ; % verify the sign!

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

F_adj = F_adj + sparse(Rrows,1,Rcoef,MESH.numNodes,1);


%% Impose zero-mean condition and solve 

A_total = [ A_adj , B ; B' , 0 ] ; 
F_total = [ F_adj ; 0 ] ; 

fprintf('\n Solve Au = f ... ');
t_solve = tic;
p                         = zeros(MESH.numNodes,1);
p_total = A_total \ F_total ; 
p(MESH.internal_dof)      = p_total(1 : end -1 ) ;
t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

%% Visualize solution 

figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',p(1:MESH.numVertices),'xystyle','interp',...
       'zdata',p(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
colormap(jet); lighting phong ; title('Solution of the adjoint problem p') 
% axis equal


