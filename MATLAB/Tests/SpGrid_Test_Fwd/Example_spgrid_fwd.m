%EXAMPLE_SPGRID_FWD solves the forward problem in the heart using the 
% stochastic collocation method

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Add path to Inverse problem folder
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')

% Add path to Sparse grid Matlab toolkit
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 

% Define which elements to use
fem = 'P1';

%% Physical domain

% meshFileName = '../../Mesh/Square/square_coarse' ;
meshFileName = '../../Mesh/Square/square_fine' ;

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

% Fill inner mesh data structure

MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_HEART_REGION ) ; 
MESH.indexInnerBoun = find( MESH.boundaries( 5 , :) == DATA.FLAG_HEART_BOUNDARY ) ;
MESH.indexInnerVert = MESH.boundaries(1, MESH.indexInnerBoun ) ;

MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;

MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 

MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% Print some information
fprintf(' * Number of Nodes           = %d \n',MESH.numNodes);
fprintf(' * Number of inner Nodes     = %d \n',MESH.numInnerNodes);




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

%% Control function

tau = 0.2 ; % smoothing level
R = 0.4 ; % radius
D = [0 0]; % displacement

w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;

w = w_target ;
% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
    H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w ) ; 
    [X,Y] = meshgrid(-1:0.02:1) ; 
    figure
    surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    shading interp ; colormap jet ; title('Target control function') ; axis equal ;
    clear H ; clear X ; clear Y ; 

%% test solveFwd and solveFwdHandler

% u = solveFwdHandler( MESH , FE_SPACE , DATA , w , [0.03 ; 3.0] ) ;
% 
% figure
% pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',u(1:MESH.numVertices),'xystyle','interp',...
%        'zdata',u(1:MESH.numVertices),'zstyle','continuous',...
%        'colorbar','on', 'mesh' , 'off'  );
% colormap(jet); lighting phong;

%% Stochastic domain 

% Here we define our stochastic space 
% Let us take for example M0 and Mi as stochastic parameters

% M0 uniformly sampled between 0.02 and 0.04
% Mi uniformly sampled between 2.0 and 3.0

N = 2 ; 
knotsM0 = @(n) knots_CC( n , 0.02 , 0.04 , 'prob' ) ; 
knotsMi = @(n) knots_CC( n , 2.0  , 4.0  , 'prob' ) ; 

level = 3 ; %level 

[ lev2knots , idxset ] = define_functions_for_rule( 'SM' , N ) ;

S = smolyak_grid(N,level,{knotsM0 , knotsMi},lev2knots , idxset); % grid

% Visualize grid
figure
plot_grid(S,[],'color','b','marker','o','MarkerFaceColor','b');

% Reduce grid
Sr = reduce_sparse_grid( S ) ;


%% Solve forward problem on the nodes

% In order to do this we define a function handler that is responsible for
% calling the appropriate function solveFwd 
% The handler is defined in the function solveFwdHandler within this file's
% directory

if ~check_if_parallel_on()
    activate_parallel() % optional argument to specify how many workers
end

U = [] ;

min_eval = 10; % minimum number of evaluations required to solve in parallel

fprintf('\nSolving laplacian on grid nodes in parallel ...\n') ;

t_evaluation = tic ; 

f = @(x) solveFwdHandler( MESH , FE_SPACE , DATA , w , x , [] , [] , [] ) ;

U = evaluate_on_sparse_grid( f , Sr , [] , [] , min_eval ) ;

t_evaluation = toc(t_evaluation);
fprintf('\nEvaluation done in %3.3f seconds.\n\n' , t_evaluation) ;


if check_if_parallel_on()
    close_parallel()
end


%% Plot some solutions 

Npl = 3;

V = randsample( 1:size(Sr.knots , 2) , Npl*Npl  ) ;


figure
for i = 1:Npl 
    for j = 1:Npl
        index = V( (i-1)*Npl + j ) ;
        subplot( Npl , Npl , (i-1)*Npl + j );
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',U(1:MESH.numVertices , index ),'xystyle','interp',...
           'zdata',U(1:MESH.numVertices , index),'zstyle','continuous',...
           'colorbar','on', 'mesh' , 'off'  );
        colormap(jet); lighting phong;
    end
end





