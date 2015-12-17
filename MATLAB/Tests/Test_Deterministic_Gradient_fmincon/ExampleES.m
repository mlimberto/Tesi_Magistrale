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

% meshFileName = '../../Mesh/Square/square_coarse' ;
% meshFileName = '../../Mesh/Square/square' ;
meshFileName = '../../Mesh/Square/square_fine' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

if (PLOT_ALL)
    figure
    pdeplot(vertices,[],elements(1:3,:))
    axis equal
end

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

%% Target solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3 : Solve the fwd problem to find zd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Evaluating target solution zd ... \n');

% Define a target control function 

tau = 0.2 ; % smoothing level
R = 0.4 ; % radius
D = [1.0 0]; % displacement

w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;

w = w_target ;
% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
if (PLOT_ALL)
    H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w_target ) ; 
    [X,Y] = meshgrid(-1:0.02:1) ; 
    figure
    surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    shading interp ; colormap jet ; title('Target control function') ; axis equal ;
    clear H ; clear X ; clear Y ; 
end

% Evaluate rhs 
F_fwd = DATA.coeffRhs * A_source_fwd * wbar ;

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

% Apply boundary conditions
fprintf('\n Apply boundary conditions ');
[A_in, F_in, u_D]   =  ApplyBC_2D(A_fwd, F_fwd, FE_SPACE, MESH, DATA);

% Impose zero-mean condition
A_total = [ A_in , B ; B' , 0 ] ; 
F_total = [ F_in ; 0 ] ; 

% Solve
fprintf('\n Solving for zd ... ');
t_solve = tic;
zd                         = zeros(MESH.numNodes,1);
zd_total = A_total \ F_total ; 
zd(MESH.internal_dof)      = zd_total(1 : end -1 ) ;
zd(MESH.Dirichlet_dof)     = u_D;t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);
clear zd_total;

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


if (PLOT_ALL)
    % control function
    H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w ) ; 
    [X,Y] = meshgrid(-1:0.02:1) ; 
    figure
    surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    shading interp ; colormap jet ; title('Control function') ; axis equal ;
    clear H ; clear X ; clear Y ; 

    % Diffusion coefficient
    [X,Y] = meshgrid( -2:0.1:2 , -2:0.1:2 ) ; 
    figure
    surf( X , Y , DATA.diffusion(X,Y) ) 
    shading interp
    colormap jet
    colorbar
    title('Diffusion coefficient')
    clear X; clear Y;
end



%% A few CG iterations

dJ_L2 = [] ;
dJ_H1 = [] ;

iterMax = 5001 ; 

% Find first search direction d0 
    % Solve fwd
    F_fwd = DATA.coeffRhs * A_source_fwd * wbar ;
    [A_in, F_in, u_D]   =  ApplyBC_2D(A_fwd, F_fwd, FE_SPACE, MESH, DATA);
    u = zeros(MESH.numNodes , 1) ;
    u_total = A_total \ [ F_in ; 0 ] ; 
    u(MESH.internal_dof) = u_total( 1 : end -1  ) ;    clear u_total; 

    % Solve adj
    F_adj = Apply_AdjBC( FE_SPACE , MESH , zd - u ) ;
    p = zeros(MESH.numNodes , 1);
    p_total = A_total' \ [ F_adj ; 0 ] ;
    p(MESH.internal_dof) = p_total(1: end -1) ;   clear p_total;

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
 

for i=1:100
    
    % Compute step length
    ACCEPTABLE = 0 ;
    
    while( ~ACCEPTABLE) 
        % Compute wnew, solve the forward and adjoint problem, evaluate new objective function and new gradient
        w_new = w + s_cg *( d  ) ;
        w_new_projected =  min( 1 , max( 0 , w_new )) ;
        w_new_projected_bar = extend_with_zero( w_new_projected , MESH) ;
            
        F_fwd = DATA.coeffRhs * A_source_fwd * w_new_projected_bar ;
        [A_in, F_in, u_D]   =  ApplyBC_2D(A_fwd, F_fwd, FE_SPACE, MESH, DATA);
        u = zeros(MESH.numNodes , 1) ;
        u_total = A_total \ [ F_in ; 0 ] ; 
        u(MESH.internal_dof) = u_total( 1 : end -1  ) ;    clear u_total; 

        F_adj = Apply_AdjBC( FE_SPACE , MESH , zd - u ) ;
        p = zeros(MESH.numNodes , 1);
        p_total = A_total' \ [ F_adj ; 0 ] ;
        p(MESH.internal_dof) = p_total(1: end -1) ;   clear p_total;

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
            
            s_cg = 1e-3;
            
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



%% Solve with constrained minimization

functionHandler = @(w) solveFwdAdjGrad( w , MESH , FE_SPACE , DATA , zd ) ;

w0 = w;
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(size(w));
ub = ones (size(w));
nonlcon = [];

options = optimoptions('fmincon','GradObj','on');

[w_opt ,J_opt ] = fmincon(functionHandler,w0,A,b,Aeq,beq,lb,ub,nonlcon,options);


%% Plot w_opt
    H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w_opt ) ; 
    [X,Y] = meshgrid(-1:0.02:1) ; 
    figure
    surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    shading interp ; colormap jet ; title('Control function') ; axis equal ;
    clear H ; clear X ; clear Y ; 

