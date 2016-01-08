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
R = 1.2 ; % radius
D = [0 6.8]; % displacement

w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;

w_target_bar = extend_with_zero( w_target , MESH ) ;

w = w_target ;
% Extend the control function to the outer boundary

wbar = extend_with_zero( w , MESH) ; 

% Visualize w
% if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',w_target_bar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis equal
    drawnow
% end

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
    
end

%% Loop 

J = [] ;

dJ_L2 = [] ;
dJ_H1 = [] ;

iterMax = 50001 ; 

errorL2 = [] ;
errorH1 = [] ;

for i=1:iterMax 
   
    % FORWARD PROBLEM 
    fprintf('\n Solving the forward problem ... ');
   
    F_fwd = DATA.coeffRhs * A_source_fwd * wbar ;
    [A_in, F_in, u_D]   =  ApplyBC_2D(A_fwd, F_fwd, FE_SPACE, MESH, DATA);
    t_solve = tic;
    u = zeros(MESH.numNodes , 1) ;
    u_total = A_total \ [ F_in ; 0 ] ; 
    u(MESH.internal_dof) = u_total( 1 : end -1  ) ;    clear u_total; 
    t_solve = toc(t_solve) ;
    fprintf('done in %3.3f s \n', t_solve); 
    
    if (PLOT_ALL)
        figure(90)
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',u(1:MESH.numVertices),'xystyle','interp',...
            'zdata',u(1:MESH.numVertices),'zstyle','continuous',...
            'colorbar','on', 'mesh' , 'off'  );
        colormap(jet); lighting phong ; title('Solution of state problem u') 
    end
    
    
    % ADJOINT PROBLEM
    fprintf('\n Solving the adjoint problem ... ');

    F_adj = Apply_AdjBC( FE_SPACE , MESH , zd - u ) ;

    t_solve = tic;
    p = zeros(MESH.numNodes , 1);
    p_total = A_total' \ [ F_adj ; 0 ] ;
    p(MESH.internal_dof) = p_total(1: end -1) ;   clear p_total;
    t_solve = toc(t_solve); fprintf('done in %3.3f s \n', t_solve);
    
    if (PLOT_ALL)
        figure(91)
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',p(1:MESH.numVertices),'xystyle','interp',...
            'zdata',p(1:MESH.numVertices),'zstyle','continuous',...
            'colorbar','on', 'mesh' , 'off'  );
        colormap(jet); lighting phong ; title('Solution of the adjoint problem p') 
    end   
    
    
    % EVALUATE OBJECTIVE FUNCTION
    J = [ J ; eval_ObjFunction(MESH , DATA , FE_SPACE , w , u , zd , -F_adj ) ] ;
    
    fprintf('\n J = %3.3f \n ',J(end) );
    
    
    % GRADIENT OF W
    fprintf('\n Solving for the gradient of w ... ');

    % Assemble source term
    F_grad = A_tBp * p ...
           + DATA.betaL2 * FE_SPACE.A_reaction_heart * wbar ...
           + DATA.betaGr * FE_SPACE.A_diffusion_heart * wbar ;
  
    t_solve = tic;
    dw = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) \ F_grad(MESH.indexInnerNodes) ;
    t_solve = toc(t_solve); fprintf('done in %3.3f s \n', t_solve); 
        
    dwbar = extend_with_zero( dw , MESH ) ; 
    
    normgradL2 = sqrt( dwbar' * FE_SPACE.A_reaction_heart * dwbar ) ;
    normgradH1 = sqrt( dwbar' * FE_SPACE.A_reaction_heart * dwbar ...
                     + dwbar' * FE_SPACE.A_diffusion_heart * dwbar ) ;
                 
    dJ_L2 = [ dJ_L2 ; normgradL2 ] ;
    dJ_H1 = [ dJ_H1 ; normgradH1 ] ;

    if (PLOT_ALL)
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',dwbar(1:MESH.numVertices),'xystyle','interp',...
       'zdata',dwbar(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    title('Control derivative Riesz element')
    end
        
    % UPDATE W
    w_new = w - DATA.gstep *( dw  ) ;
    % Apply projection step
    w =  min( 1 , max( 0 , w_new )) ;
    
    wbar = extend_with_zero( w , MESH) ;
    
    if (PLOT_ALL)
    figure(92)
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
    'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis equal
    drawnow
    end    
   
    % Save a snapshot every once in a while ... 
    if mod( i , 500 ) == 1 
        figure
        pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',wbar(1:MESH.numVertices),'xystyle','interp',...
            'zdata',wbar(1:MESH.numVertices),'zstyle','continuous',...
            'colorbar','on', 'mesh' , 'off' );
        colormap(jet);
        lighting phong
        view([0 90])
        axis equal
        drawnow 
    end

%     figure(66)
%     b = MESH.boundaries(1:2 , find( MESH.boundaries(5,: ) ~= 3  ) ) ;
%     plot3(MESH.vertices(1,b(1,:)) , MESH.vertices(2,b(1,:)) ,zd(b(1,:)) , 'Linewidth',2)
%     hold on
%     plot3(MESH.vertices(1,b(1,:)) , MESH.vertices(2,b(1,:)) ,u(b(1,:)) , 'Linewidth',2)
%     legend('zd','u')
%     drawnow
%     hold off
    
% Evaluate the error

errorL2 = [ errorL2 ; sqrt(productL2Heart( w - w_target , w - w_target , MESH , FE_SPACE ) ) ] ;
errorH1 = [ errorH1 ; sqrt(productH1Heart( w - w_target , w - w_target , MESH , FE_SPACE ) ) ] ;

end

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

if (PLOT_ALL)
figure
outer_oundary_index = MESH.boundaries(1:2 , find( MESH.boundaries(5,: ) ~= 3  ) ) ;
plot3(MESH.vertices(1,outer_oundary_index(1,:)) , MESH.vertices(2,outer_oundary_index(1,:)) ,zd(outer_oundary_index(1,:)) , 'Linewidth',2)
hold on
plot3(MESH.vertices(1,outer_oundary_index(1,:)) , MESH.vertices(2,outer_oundary_index(1,:)) ,u(outer_oundary_index(1,:)) , 'Linewidth',2)
legend('zd','u')
figure
plot3(MESH.vertices(1,outer_oundary_index(1,:)) , MESH.vertices(2,outer_oundary_index(1,:)) ,abs(zd(outer_oundary_index(1,:)) - u(outer_oundary_index(1,:)) ) , 'Linewidth',2)
legend('u-zd at the boundary')
end

