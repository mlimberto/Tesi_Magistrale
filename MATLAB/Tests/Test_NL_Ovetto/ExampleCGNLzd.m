%EXAMPLE solves the deterministic control problem.

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')

% Plotting settings
PLOT_ALL = 0 ; % 1 or 0 value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1 : Mesh and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define which elements to use
fem = 'P1';

%% Import mesh

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

%% Set parameters from file
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

% Initialize empty matrices in the FE_SPACE structure

FE_SPACE.A_diffusion_heart = [] ;
FE_SPACE.A_reaction_heart = [] ;
FE_SPACE.A_diffusion_total = [] ;

%% Matrix assembly

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

DATA_temp = DATA ;
DATA_temp.diffusion = @(x,y,t,param) 1 + 0*x.*y ;

% Assemble rhs matrix for forward problem
fprintf('\n Assembling source term matrix ... ');
t_assembly_source = tic;
A_source_fwd          =  Assembler_2D(MESH, DATA_temp, FE_SPACE , 'diffusion' , [] , [] , DATA.FLAG_HEART_REGION );
t_assembly_source = toc(t_assembly_source);
fprintf('done in %3.3f s\n', t_assembly_source);

clear DATA_temp ;
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

% D = [-3.9 0]; % displacement
% 
% w_target2 = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
% w_target2 = 1 - smoothLS(w_target2 , tau) ;

w_target_bar = extend_with_zero( w_target , MESH )  ... 
%                 + extend_with_zero( w_target2 , MESH )  ... 
                ;

w = w_target ...
%     + w_target2 ...
    ;

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

if (PLOT_ALL)
% Visualize rhs term 
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',F_fwd(1:MESH.numVertices),'xystyle','interp',...
        'zdata',F_fwd(1:MESH.numVertices),'zstyle','continuous',...
        'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    title('Source term')
end


F_fwd = DATA.coeffRhs * A_source_fwd * extend_with_zero(w,MESH) ;
zd = solveFwdNL( MESH , FE_SPACE , DATA , w , [] , B , F_fwd ) ;


if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
       'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off'  );
    colormap(jet); lighting phong; title('Target solution z_d')
%     axis equal
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4 : Iterative method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize a control function

% clear previous data
clear w ; clear wbar ;

tau = 0.2 ; % smoothing level

% Since we use a level set we define radius and center coordinates
R = 0.4 ; % radius
D = [0 0]; % displacement

% w0 = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
% w0 = 1 - smoothLS(w0 , tau) ;

w0 = zeros(MESH.numInnerNodes , 1);
% w0 = ones(MESH.numInnerNodes , 1);

w = w0 ;
wbar = extend_with_zero( w , MESH) ; 


%% Loop initialization

dJ_L2 = [] ;
dJ_H1 = [] ;

error_L2 = [];
error_H1 = [];

iterMax = 5001 ; 

% Find first search direction d0 
    % Solve fwd
    F_fwd = DATA.coeffRhs * A_source_fwd * wbar ;
    
    u = solveFwd( MESH , FE_SPACE , DATA , w , A_fwd , FE_SPACE.B , F_fwd ) ;
    
    [ p , F_adj ] = solveAdj( MESH , FE_SPACE , DATA , w , u , zd, A_fwd' , FE_SPACE.B  ) ;
    
    % Find gradient of J and set first descent direction as -grad(J)
    F_grad = A_tBp * p ...
          + DATA.betaL2 * FE_SPACE.A_reaction_heart * wbar ...
          + DATA.betaGr * FE_SPACE.A_diffusion_heart * wbar ;
  
    dw =  A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ F_grad(MESH.indexInnerNodes)  ;
    
    d = - dw ;
    dbar = extend_with_zero( d , MESH ) ; 
        
   
% Set an initial step length
s_cg = 1e-2; 
   
% Set Wolfe coefficients 
sigma1_cg = 1e-7 ; 
sigma2_cg = 0.4;
 
% Compute first L2 norm of gradient (indeed it's the L2 norm squared)
normgradL2 = productL2Heart( dw , dw , MESH , FE_SPACE ) ;
normgradL2_old = 0 ;

normgradH1 = productH1Heart( dw , dw , MESH , FE_SPACE ) ;
normgradH1_old = 0;

% Evaluate first objective function 
J_old = eval_ObjFunction(MESH , DATA , FE_SPACE , w , u , zd , -F_adj ) ;
J = [ J_old ] ;

% Initialize d_old
d_old = d ;

%% Loop 

for i=1:2000
    
    % Compute step length
    ACCEPTABLE = 0 ;
    
    while( ~ACCEPTABLE) 
        % Compute wnew, solve the forward and adjoint problem, evaluate new objective function and new gradient
        w_new = w + s_cg *( d  ) ;
        w_new_projected =  min( 1 , max( 0 , w_new )) ;
        w_new_projected_bar = extend_with_zero( w_new_projected , MESH) ;
            
        F_fwd = DATA.coeffRhs * A_source_fwd * w_new_projected_bar ;
        u = solveFwd( MESH , FE_SPACE , DATA , w_new_projected , A_fwd , FE_SPACE.B , F_fwd ) ;
    
        [ p , F_adj ] = solveAdj( MESH , FE_SPACE , DATA , w_new_projected , u , zd, A_fwd' , FE_SPACE.B  ) ;


        F_grad = A_tBp * p ...
               + DATA.betaL2 * FE_SPACE.A_reaction_heart * wbar ...
               + DATA.betaGr * FE_SPACE.A_diffusion_heart * wbar ;
  
        dw_new = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ F_grad(MESH.indexInnerNodes) ;
        
        % Evaluate new J 
        J_new = eval_ObjFunction(MESH , DATA , FE_SPACE , w_new_projected , u , zd , -F_adj ) ;
        
        % Check Armijo's condition (Wolfe condition 1)
        W1 = J_new <= J_old + sigma1_cg * s_cg * productH1Heart( dw , d , MESH , FE_SPACE ) ;
        
        % Check curvature condition (Wolfe condition 2)
        W2 = productH1Heart( dw_new , d , MESH , FE_SPACE ) >= ...
             sigma2_cg * productH1Heart( dw , d , MESH , FE_SPACE ) ;
         
        W2 = 1 ; 
        
        if ( s_cg < 1e-7 ) 
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
            
            s_cg = 1e-2;
            
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
    
    error_L2 = [ error_L2 ; sqrt( productL2Heart(w - w_target , w - w_target , MESH , FE_SPACE ) )  ] ;
    error_H1 = [ error_H1 ; sqrt( productH1Heart(w - w_target , w - w_target , MESH , FE_SPACE ) )  ] ;

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

%% Final visualizations

% Plot solution
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off' );
colormap(jet);
lighting phong
view([0 90])
axis equal
drawnow

% Display objective function 
% if (PLOT_ALL)
figure
loglog( J , 'LineWidth' , 2 ) 
hold on
loglog( dJ_L2, 'LineWidth' , 2  ) ; 
loglog( dJ_H1 , 'LineWidth' , 2 ) ;
legend('j','l2','h1')
grid on 
% end

% Plot error
figure
subplot(1,2,1)
loglog( error_L2, 'LineWidth' , 2  ) ; 
legend('L2 error')
grid on 
% figure
subplot(1,2,2)
loglog( error_H1 , 'LineWidth' , 2 ) ;
legend('H1 error')
grid on 

if (PLOT_ALL)
% This plot has to be fixed
figure
b = MESH.boundaries(1:2 , find( MESH.boundaries(5,: ) == 13  ) ) 
plot3(MESH.vertices(1,b(1,:)) , MESH.vertices(2,b(1,:)) ,zd(b(1,:)) , 'Linewidth',2)
hold on
plot3(MESH.vertices(1,b(1,:)) , MESH.vertices(2,b(1,:)) ,u(b(1,:)) , 'Linewidth',2)
legend('zd','u')

figure
plot3(MESH.vertices(1,b(1,:)) , MESH.vertices(2,b(1,:)) ,abs(zd(b(1,:)) - u(b(1,:)) ) , 'Linewidth',2)
legend('u-zd at the boundary')
end


if (PLOT_ALL)
% Plot solution vs target solution
figure
subplot(1,2,1)
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','off', 'mesh' , 'off' );
colormap(jet);
lighting phong
view([0 90])
axis equal
subplot(1,2,2)
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','off', 'mesh' , 'off' );
colormap(jet);
lighting phong
view([0 90])
axis equal
end

if (PLOT_ALL)
% Plot heart and torso in a nice colour
figure
a = ones( size( wbar ) ) ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',a(1:MESH.numVertices),'xystyle','flat',...
    'colorbar','off', 'mesh' , 'off' );
hold on 
a = zeros( size( wbar ) ) ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,find( MESH.elements(end,:) == 2 ) ),'xydata',a(1:MESH.numVertices),'xystyle','flat',...
    'colorbar','off', 'mesh' , 'off' );
view([0 90])
axis equal
end
