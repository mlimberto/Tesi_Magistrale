%% Add paths

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')
addpath('../../Inverse_Problem/')
addpath('../../Inverse_Problem/Examples/')
addpath( genpath( '../../sparse-grids-matlab-kit/' ) ) 
addpath( genpath( '../../RF' ) )

%% Create mesh

addpath('../../Mesh/Ovetto_Circondato/') ; % add path to mesh geometry specification function
[vertices,boundaries,elements] = initmesh('ovetto_circondato_sens' ,'Jiggle','minimum','Hgrad',1.9,'Hmax',1000) ;
% Refine only the heart subdomain
[vertices,boundaries,elements] = refinemesh('ovetto_circondato_sens' , vertices , boundaries, elements , [1 ]) ;                                    
% Refine both the heart and the torso                                        
[vertices,boundaries,elements] = refinemesh('ovetto_circondato_sens' , vertices , boundaries, elements , [1 2 ]) ;                                         
  
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
% BUG evidente nel loop qui sopra se si sceglie ad esempio 6,5,13 !!!


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

% Set finite elements to use
fem = 'P1' ;

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

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf2D(MESH.vertices, MESH.elements);   

% Compute quadrature nodes and weights on the reference element
quad_order                  = 4; 
[quad_nodes, quad_weights]  = dunavant_quad(quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis2D('P1', quad_nodes(1,:), quad_nodes(2,:));

% Fill inner mesh data structure

MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_HEART_REGION ) ; 
MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;
MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 
MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% Print some information
fprintf(' * Number of Nodes           = %d \n',MESH.numNodes);
fprintf(' * Number of inner Nodes     = %d \n',MESH.numInnerNodes);

%% Set random field covariance structure

corr.name = 'gauss';
corr.c0 = [ 2 1] ; 

% corr.sigma = cos(pi*MESH.innerNodes(1,:)').*sin(2*pi*MESH.innerNodes(2,:)')+1.5;
corr.sigma = repmat(1.5 , size(MESH.innerNodes,2) ,1 );

% Generate RF expansion 

N = 20 ;

[F,KL] = randomfield(corr,MESH.innerNodes','trunc', N);

%% Get a realization and plot it

N_plots = 2 ; 

NS = N_plots^2 ; 

F = repmat(KL.mean,1,NS) + KL.bases*diag(KL.sv)*randn(N,NS) ;

figure

for i=1:NS
    F_bar = extend_with_zero( F(:,i) , MESH ); 
    subplot(N_plots,N_plots,i);
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',F_bar(1:MESH.numVertices),'xystyle','interp',...
        'zdata',F_bar(1:MESH.numVertices),'zstyle','continuous',...
        'colorbar','on', 'mesh' , 'off'  );
    colormap(jet); lighting phong; title('Realization of the random field')
    view([0 90])
    axis equal
end










