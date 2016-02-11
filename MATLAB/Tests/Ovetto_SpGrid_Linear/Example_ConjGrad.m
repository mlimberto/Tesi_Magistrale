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
D = [0 +6.8]; % displacement

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

zd = solveFwd( MESH , FE_SPACE , DATA , w  ) ; 

% Visualize zd
% if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
       'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
    colormap(jet); lighting phong; title('Target solution z_d')
%     axis equal
% end

%% Initialize a control function

% clear previous data
clear w ; clear wbar ;

fprintf('\n Initializing control function ... \n');

w0 = zeros(MESH.numInnerNodes , 1);
% w0 = ones(MESH.numInnerNodes , 1);

w = w0 ;
wbar = extend_with_zero( w , MESH) ; 

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

% Assemble rhs matrix for forward problem
temp = DATA.diffusion ;
DATA.diffusion = @(x,y,t,param) 1 + 0*x.*y ;
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
FE_SPACE.A_diffusion_heart          =  Assembler_2D(MESH, DATA, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);
DATA.diffusion = temp ;

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
FE_SPACE.A_reaction_heart = A_grad_L2 ;
% diffusion part ( i.e. H_0^1 product has already been assembled and it's called A_source_fwd)
A_grad_H01 = FE_SPACE.A_diffusion_heart ; 
A_grad = A_grad_L2 + A_grad_H01 ; 
t_assembly_grad = toc(t_assembly_grad);  fprintf('done in %3.3f s\n', t_assembly_grad);


%% Solve forward and adjoint problem on the nodes

if ~check_if_parallel_on()
    activate_parallel() % optional argument to specify how many workers
end

min_eval = 10; % minimum number of evaluations required to solve in parallel

% To close parallel pool then type the following
% if check_if_parallel_on()
%     close_parallel()
% end

%% Loop initialization

dJ_L2 = [] ;
dJ_H1 = [] ;

error_L2 = [];
error_H1 = [];

iterMax = 5001 ; 

% Find first search direction d0

    % Solve forward and adjoint problem    
    fprintf('\n Solving fwd and adjoint ... ');
    t_solve = tic;
    f = @(x) solveFwdAdjHandler( MESH , FE_SPACE , DATA , w , zd , x  ) ;
%     [UP] = evaluate_on_sparse_grid( f , Sr , [] , []  ) ; % Serial
    [UP] = evaluate_on_sparse_grid( f , Sr , [] , [] , min_eval ) ; % Parallel
    t_solve = toc(t_solve); fprintf('done in %3.3f s \n', t_solve); 

    F_grad = UP(1:MESH.numNodes , :) ;
    U = UP( (MESH.numNodes+1) : (2*MESH.numNodes) , : ) ;
    P = UP( (2*MESH.numNodes+1) : end -1 , : ) ;
    
    % EVALUATE OBJECTIVE FUNCTION
    J_old =  UP(end,:) * Sr.weights'  ;
    J = [ J_old ] ;
    fprintf('\n J = %3.3f \n ',J(end) );
    
    % Compute gradient    
    
    E_F_grad = F_grad * Sr.weights' ;
    dw = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ E_F_grad(MESH.indexInnerNodes) ;
    
    d = - dw ;
    dbar = extend_with_zero( d , MESH ) ; 
        
   
% Set an initial step length
s_cg_init = 1e-3;
s_cg = s_cg_init; 
   
% Set Wolfe coefficients 
sigma1_cg = 1e-7 ; 
sigma2_cg = 0.4;
 
% Compute first L2 norm of gradient (indeed it's the L2 norm squared)
normgradL2 = productL2Heart( dw , dw , MESH , FE_SPACE ) ;
normgradL2_old = 0 ;

normgradH1 = productH1Heart( dw , dw , MESH , FE_SPACE ) ;
normgradH1_old = 0;

% Initialize d_old
d_old = d ;



%% Loop 

for i=1:iterMax
    
    % Compute step length
    ACCEPTABLE = 0 ;
    
    while( ~ACCEPTABLE) 
        % Compute wnew, solve the forward and adjoint problem, evaluate new objective function and new gradient
        w_new = w + s_cg *( d  ) ;
        w_new_projected =  min( 1 , max( 0 , w_new )) ;
        w_new_projected_bar = extend_with_zero( w_new_projected , MESH) ;
            
        % Solve forward and adjoint problem    
        fprintf('\n Solving fwd and adjoint ... ');
        t_solve = tic;
        f = @(x) solveFwdAdjHandler( MESH , FE_SPACE , DATA , w_new_projected , zd , x  ) ;
%     [UP] = evaluate_on_sparse_grid( f , Sr , [] , []  ) ; % Serial
        [UP] = evaluate_on_sparse_grid( f , Sr , [] , [] , min_eval ) ; % Parallel
        t_solve = toc(t_solve); fprintf('done in %3.3f s \n', t_solve); 

        F_grad = UP(1:MESH.numNodes , :) ;
        U = UP( (MESH.numNodes+1) : (2*MESH.numNodes) , : ) ;
        P = UP( (2*MESH.numNodes+1) : end -1 , : ) ;
    
        % EVALUATE OBJECTIVE FUNCTION
        J_new =   UP(end,:) * Sr.weights'  ;
%         fprintf('\n J = %3.3f \n ',J(end) );
    
        % Compute gradient    
        E_F_grad = F_grad * Sr.weights' ;
        dw_new = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ E_F_grad(MESH.indexInnerNodes) ;
  
        % Check Armijo's condition (Wolfe condition 1)
        W1 = J_new <= J_old + sigma1_cg * s_cg * productH1Heart( dw , d , MESH , FE_SPACE ) ;
        
        % Check curvature condition (Wolfe condition 2)
        W2 = productH1Heart( dw_new , d , MESH , FE_SPACE ) >= ...
             sigma2_cg * productH1Heart( dw , d , MESH , FE_SPACE ) ;
         
        W2 = 1 ; 
        
        if ( s_cg < 1e-9 ) 
           W1 = 1 ; 
           d = -dw_new ; 
           d_old = -dw_new ;
           disp 'Reinitializing direction'
        end
        
        % If conditions are verified update the solutions
        % Otherwise, update the gradient step
        if ( W1 && W2 ) 
            ACCEPTABLE = 1 ;
            
            w = min( 1 , max( 0 , w_new )) ;
            wbar = extend_with_zero( w , MESH) ;
            J = [ J ; J_new ] ;
            J_old = J_new ;
            dw_old = dw ;
            dw = dw_new ;
            dwbar = extend_with_zero( dw , MESH ) ;
            
            s_cg = s_cg_init; 

            fprintf('\n J = %3.3f \n ',J(end) );
        else

            if ~W1
                s_cg = s_cg / 1.2 ;
                disp 'Decreasing step'
            end

            
        end
        
    end
                
    % Evaluate L2 and H1 norm of new gradient
    normgradL2_old = normgradL2 ;
    normgradL2 = productL2Heart( dw , dw , MESH , FE_SPACE ) ;
    normgradH1_old = normgradH1 ;
    normgradH1 = productH1Heart( dw , dw , MESH , FE_SPACE ) ;
    
    
    dJ_L2 = [ dJ_L2 ; sqrt( normgradL2 )  ] ;
    dJ_H1 = [ dJ_H1 ; sqrt( normgradH1 )  ] ;
    
%     error_L2 = [ error_L2 ; sqrt( productL2Heart(w - w_target , w - w_target , MESH , FE_SPACE ) )  ] ;
%     error_H1 = [ error_H1 ; sqrt( productH1Heart(w - w_target , w - w_target , MESH , FE_SPACE ) )  ] ;
        
    % Determine a scalar beta
    
    % Fletcher-Reeves rule
    beta_cg_FR = normgradH1 / normgradH1_old ; 
    
    % Polak-Ribiere rule
    beta_cg_PR = productH1Heart( dw , dw - dw_old , MESH , FE_SPACE ) / normgradH1 ;  
    
    % Hestenes-Stiefel rule
    if ( i == 1 )
    beta_cg_HS = beta_cg_PR ;       
    else
    beta_cg_HS = productH1Heart( dw , dw - dw_old , MESH , FE_SPACE ) ...
               / productH1Heart( d_old , dw - dw_old , MESH , FE_SPACE ) ;    
    end
    
    % Find the new descent direction 
    d_old = d ;
    
    d =  -dw + beta_cg_FR * d ; 
%     d =  -dw + beta_cg_PR * d ; 
%     d = -dw + beta_cg_HS * d ;

    % Check termination conditions
        
        % Check if d is a descent direction
        d_dot_dw = productH1Heart( d , dw , MESH , FE_SPACE ) / sqrt(normgradH1) / sqrt( productH1Heart( d , d , MESH , FE_SPACE ) ) ;
        fprintf('\n cos(angle) = %3.3f \n ',d_dot_dw );
end

%% Visualisations

        % Optimal control
        figure
        pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
            'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
            'colorbar','on', 'mesh' , 'off' );
        colormap(jet);
        lighting phong
        view([0 90])
        axis equal
        drawnow 
        
        
        figure
        loglog(J , 'LineWidth',2 ) ; 
        grid on ;










