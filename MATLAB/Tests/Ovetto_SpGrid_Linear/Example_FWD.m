% This file solves the forward problem on an oval mesh using 
% stochastic collocation for discretization in the stochastic space 

clc
clear all
% Add path to FE solver, Inverse problem folder and SP grids toolkit
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 

% Define which elements to use
fem = 'P1';

% Plotting settings
PLOT_ALL = 0 ;

%% Setting up the spatial part

addpath('../../Mesh/Ovetto_Circondato/') ; % add path to mesh geometry specification function

[vertices,boundaries,elements] = initmesh('ovetto_circondato_sens' , ...
                                          'Jiggle','minimum','Hgrad',1.9,'Hmax',1000) ;

% Refine only the heart subdomain
[vertices,boundaries,elements] = refinemesh('ovetto_circondato_sens' , vertices , ... 
                                            boundaries, elements , [1 ]) ;                                    

% Refine both the heart and the torso                                        
[vertices,boundaries,elements] = refinemesh('ovetto_circondato_sens' , vertices , ... 
                                            boundaries, elements , [1 2 ]) ;                                         
  
% Reinitialize the boundary indexes
for i = 1:4 
    boundaries( 5 , find( boundaries(5,:) == i)  ) = 3 ;  
end
for i = 5:12
    boundaries( 5 , find( boundaries(5,:) == i)  ) = 5 ;
end
for i = 13:26
    boundaries( 5 , find( boundaries(5,:) == i)  ) = 13 ;
end                                        

% ATTENZIONE, i valori 3,5,13 non sono scelti completamente a caso,
% BUG evidente se si sceglie ad esempio 6,5,13 !!!

if (PLOT_ALL)
    % Plot mesh
    figure
    pdeplot(vertices,[],elements(1:3,:))
    axis equal    
end

% Set parameters from file
data_file = 'DirichletNeumann_data' ;

DATA   = read_DataFile(data_file);
DATA.param = [] ;
t = [] ;

% Add flags to identify mesh elements (in accordance with .geo file)
DATA.FLAG_HEART_REGION = 1 ;
DATA.FLAG_TORSO_REGION = 2 ;

DATA.FLAG_HEART_OUTER_BOUNDARY = 3 ;
DATA.FLAG_HEART_INNER_BOUNDARY = 5 ;
DATA.FLAG_TORSO_BOUNDARY_NEU = 13 ;

% Fill MESH data structure
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
% MESH.indexInnerBoun = find( MESH.boundaries( 5 , :) == DATA.FLAG_HEART_BOUNDARY ) ;
% MESH.indexInnerVert = MESH.boundaries(1, MESH.indexInnerBoun ) ;

MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;

MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 

MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% Print some information
fprintf(' * Number of Nodes           = %d \n',MESH.numNodes);
fprintf(' * Number of inner Nodes     = %d \n',MESH.numInnerNodes);


% Plot some useful stuff
if (PLOT_ALL)
   figure
   subplot(1,2,1)
   pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem))
   axis equal   
   subplot(1,2,2)
   pdeplot(MESH.vertices,[],MESH.elements(1:3,find(MESH.elements(end,:) == 2 ) ) )
   axis equal
end

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

% Initialize empty matrices in the FE_SPACE structure

FE_SPACE.A_diffusion_heart = [] ;
FE_SPACE.A_reaction_heart = [] ;
FE_SPACE.A_diffusion_total = [] ;


%% Target solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3 : Solve the fwd problem to find zd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Evaluating target solution zd ... \n');

% Define a target control function 

tau = 0.2 ; % smoothing level
R = 1.2 ; % radius
D = [-3.9 0]; % displacement

w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;

% D = [0 6.8]; % displacement
% 
% w_target2 = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
% w_target2 = 1 - smoothLS(w_target2 , tau) ;

w_target_bar = extend_with_zero( w_target , MESH )  ...
%                 + extend_with_zero( w_target2 , MESH )  ... 
                ;

w = w_target ...
%     + w_target2 ...
    ;

% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
% if (PLOT_ALL)
    fig_target = figure ;
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis equal
    drawnow
% end


%% Setting up the Stochastic part 

% Here we define our stochastic space 
% Let us take for example M0 and Mi as stochastic parameters

% M0 uniformly sampled between 2.30 and 2.50
% Mi uniformly sampled between 2.9 and 3.1

N = 2 ; 
knotsM0 = @(n) knots_CC( n , 2.30 , 2.50 , 'prob' ) ; 
knotsMi = @(n) knots_CC( n , 2.90  , 3.10  , 'prob' ) ; 

level = 3 ; %level 

[ lev2knots , idxset ] = define_functions_for_rule( 'SM' , N ) ;

S = smolyak_grid(N,level,{knotsM0 , knotsMi},lev2knots , idxset); % grid

% Visualize grid
figure
plot_grid(S,[],'color','b','marker','o','MarkerFaceColor','b');

% Reduce grid
Sr = reduce_sparse_grid( S ) ;

%% Pre-assemble stuff

% Build vector for zero-mean condition
fprintf('\n Assembling zero-mean vector ... ');
t_assembly_zeroMean = tic;
A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
B = A_react * ones( MESH.numNodes , 1 ) ; 
t_assembly_zeroMean = toc(t_assembly_zeroMean);
fprintf('done in %3.3f s\n', t_assembly_zeroMean);
clear A_react ;

FE_SPACE.B = B ; 

%% Solve forward problem on the nodes

% In order to do this we define a function handler that is responsible for
% calling the appropriate function solveFwd 
% The handler is defined in the function solveFwdHandler within this file's
% directory

% if ~check_if_parallel_on()
%     activate_parallel() % optional argument to specify how many workers
% end

U = [] ;

min_eval = 10; % minimum number of evaluations required to solve in parallel

fprintf('\nSolving forward problem on grid nodes in parallel ...\n') ;

t_evaluation = tic ; 

f = @(x) solveFwdHandler( MESH , FE_SPACE , DATA , w , x , [] , FE_SPACE.B , [] ) ;

U = evaluate_on_sparse_grid( f , Sr , [] , []  ) ; % Serial
% U = evaluate_on_sparse_grid( f , Sr , [] , [] , min_eval ) ; % Parallel


t_evaluation = toc(t_evaluation);
fprintf('\nEvaluation done in %3.3f seconds.\n\n' , t_evaluation) ;


% if check_if_parallel_on()
%     close_parallel()
% end


%% Plot some solutions 

Npl = 3;

V = randsample( 1:size(Sr.knots , 2) , Npl*Npl  ) ;

% Set color limits
climit = [ min(min( U(: , V ) )) , max(max( U(:,V) ) ) ] ;

figure
for i = 1:Npl 
    for j = 1:Npl
        index = V( (i-1)*Npl + j ) ;
        subplot( Npl , Npl , (i-1)*Npl + j );
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',U(1:MESH.numVertices , index ),'xystyle','interp',...
           'zdata',U(1:MESH.numVertices , index),'zstyle','continuous',...
           'colorbar','on', 'mesh' , 'off'  );
        colormap(jet); lighting phong;
        caxis(climit) ;
        title_str=  sprintf('Potential for M0 = %1.3f , Mi = %3.1f' , Sr.knots(1 ,index) , Sr.knots(2 , index ) ) ;
        title(title_str);
    end 
end
