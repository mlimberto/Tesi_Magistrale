%% Example ovetto with random field for M0. No uncertainty on the other parameters.

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 
addpath( genpath( '../../RF' ) )

PLOT_ALL = 0 ;

% define which elements to use

fem = 'P1'; 

%% Import mesh and DataFile and generate FE information

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
MESH.indexOuterElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_TORSO_REGION ) ; 

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

% Initialize empty matrices in the FE_SPACE structure

FE_SPACE.A_diffusion_heart = [] ;
FE_SPACE.A_reaction_heart = [] ;
FE_SPACE.A_diffusion_total = [] ;


%% Initial Matrix assembly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 : Matrix assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assemble stiffness matrix for forward problem
fprintf('\n Assembling state matrix ... ');
t_assembly_fwd = tic;
A_fwd              =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion');
t_assembly_fwd = toc(t_assembly_fwd);
fprintf('done in %3.3f s\n', t_assembly_fwd);

FE_SPACE.A_diffusion_total = A_fwd ;

% Fix diffusion coefficient!!! 

DATA.diffusion = @(x,y,t,param) 1 + 0*x.*y ;

% Assemble rhs matrix for forward problem
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
A_source_fwd          =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);

FE_SPACE.A_diffusion_heart = A_source_fwd ;

% Build vector for zero-mean condition
fprintf('\n Assembling zero-mean vector ... ');
t_assembly_zeroMean = tic;
A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' ) ; 
B = A_react * ones( MESH.numNodes , 1 ) ; 
t_assembly_zeroMean = toc(t_assembly_zeroMean);
fprintf('done in %3.3f s\n', t_assembly_zeroMean);
clear A_react ;

FE_SPACE.B = B ; 

% Assemble lhs matrix for gradient computation
fprintf('\n Assembling lhs matrix for gradient computation ... ');
t_assembly_grad = tic;
% reaction part ( i.e. L2 product )
A_grad_L2 = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' , [] , [] , DATA.FLAG_HEART_REGION ) ; 
% diffusion part ( i.e. H_0^1 product has already been assembled and it's called A_source_fwd)
A_grad_H01 = A_source_fwd ; 
A_grad = A_grad_L2 + A_grad_H01 ; 
t_assembly_grad = toc(t_assembly_grad);  fprintf('done in %3.3f s\n', t_assembly_grad);
% clear A_grad_H01 ; clear A_grad_L2 ;

FE_SPACE.A_reaction_heart = A_grad_L2 ;

% Assemble rhs gradient matrix
    % We wish to assemble the matrix evaluating the source term B' * p 
    % We can obtain it again from the A_source_fwd matrix
A_tBp = - DATA.coeffRhs * A_source_fwd ; 

fprintf('\nTotal assembling time : %3.3f s\n', t_assembly_fwd + t_assembly_source + t_assembly_zeroMean + t_assembly_grad);


%% Define a target control function wbar

fprintf('\n Evaluating target solution zd ... \n');

% Define a target control function 

tau = 0.2 ; % smoothing level
R = 1.2 ; % radius
D = [-3.9 0]; % displacement

w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;

w_target_bar = extend_with_zero( w_target , MESH )  ;

w = w_target ;

% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
% if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis equal
    drawnow
% end

%% TEMPORARY : COMPUTE RANDOM SOLUTION WITHOUT RANDOM FIELD

% Change parameters for target solution
DATA_Zd = DATA ;
DATA_Zd.M0 = 2.39 ;
DATA_Zd.Mi = 3.0 ;
DATA_Zd.diffusion = @(x,y,t,param)( DATA_Zd.M0 + (DATA_Zd.Mi + DATA_Zd.Me - DATA_Zd.M0)*( 1 - smoothLS( DATA_Zd.heartLS(x,y) , DATA_Zd.tauDiff) ) + 0.*x.*y);
DATA_Zd.coeffRhs = -1. * (DATA_Zd.vTr_i - DATA_Zd.vTr_e )*DATA_Zd.Mi ;

% E' necessario ridefinire il termine diffusion perch? ? stato
% precedentemente modificato. Non molto leggibile come cosa.

% Solve
fprintf('\n Solving for zd ... ');
t_solve = tic ;
zd_NO_RF = solveFwdNL( MESH , FE_SPACE , DATA_Zd , w ) ; 
t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

% if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd_NO_RF(1:MESH.numVertices),'xystyle','interp',...
       'zdata',zd_NO_RF(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
    colormap(jet); lighting phong; 
%     title('Target solution z_d')
%     axis equal
% end



%% Generate random field for M0

% Set RF covariance structure
corr.name = 'gauss';
corr.c0 = [ 30 30] ; 
corr.sigma = 0.4* ones(size(MESH.nodes,2) ,1 );

% Generate RF expansion 

N = 10 ; % number of KL bases

[F,KL] = randomfield(corr,MESH.nodes','trunc', N , 'mean' , 2.39);

% Get a realization and plot it 
NS = 1 ; % number of samples

sample = randn(N,NS) ;

F_RF = repmat(KL.mean,1,NS) + KL.bases*diag(KL.sv)*sample ;

% if (PLOT_ALL)
   figure
   pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexOuterElem),'xydata',F_RF(1:MESH.numVertices),'xystyle','interp',...
        'zdata',F_RF(1:MESH.numVertices),'zstyle','continuous',...
        'colorbar','on', 'mesh' , 'off'  );
   colormap(jet); lighting phong; title('Realization of the random field')
   view([0 90])
   axis equal
% end

% UTILE INSERIRE INDEX OUTER ELEMENTS??

%% Assembla matrici di diffusione e calcola zd

% Change parameters for target solution
DATA_Zd = DATA ;
DATA_Zd.M0 = 2.39 ;
DATA_Zd.Mi = 3.0;
DATA_Zd.diffusion = @(x,y,t,param)( DATA_Zd.M0 + (DATA_Zd.Mi + DATA_Zd.Me - DATA_Zd.M0)*( 1 - smoothLS( DATA_Zd.heartLS(x,y) , DATA_Zd.tauDiff) ) + 0.*x.*y);
DATA_Zd.coeffRhs = -1. * (DATA_Zd.vTr_i - DATA_Zd.vTr_e )*DATA_Zd.Mi ;

% Parte sul torso che dipende dal RF
% Example call : A = Assembler_2D_RF(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , FLAG_HEART_REGION, [] , F )
A_rf = NonLinear_Assembler_2D_RF( MESH , DATA_Zd , FE_SPACE  , 'diffusion' , [] , [] , DATA.FLAG_TORSO_REGION  , [] , F_RF ) ;

% Parte sul cuore (in questo caso indipendente dal RF)
A_internal = Assembler_2D(MESH, DATA_Zd, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );

% Parte sul cuore che dipende dal controllo w
A_wdep  =  NonLinear_Assembler_2D(MESH, DATA_Zd, FE_SPACE , 'diffusion' , ...
                                  [] , [] , DATA.FLAG_HEART_REGION ,[] ,w );


% Calcola zd
fprintf('\n Solving for zd ... ');
t_solve = tic ;
zd = solveFwdNL( MESH , FE_SPACE , DATA_Zd , w , A_internal + A_rf + A_wdep , B  ) ; 
t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

% if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
       'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
    colormap(jet); lighting phong; 
%     title('Target solution z_d')
%     axis equal
% end



%% 








