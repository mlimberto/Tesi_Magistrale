% Define a control function only on the internal domain

clc
clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Define which elements to use
fem = 'P2';

%% Import mesh

meshFileName = '../../Mesh/Square/square_coarse' ;
% meshFileName = '../../Mesh/Square/square' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

% Add flags to identify mesh elements (in accordance with .geo file)

FLAG_HEART_REGION = 10 ;
FLAG_TORSO_REGION = 11 ;

FLAG_TORSO_BOUNDARY_DIRI = 1 ; 
FLAG_TORSO_BOUNDARY_NEU = 2 ;
FLAG_HEART_BOUNDARY = 3 ;

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


%% Fill complete MESH data structure
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


%% Fill inner mesh data structure

MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== FLAG_HEART_REGION ) ; 
MESH.indexInnerBoun = find( MESH.boundaries( 5 , :) == FLAG_HEART_BOUNDARY ) ;
MESH.indexInnerVert = MESH.boundaries(1, MESH.indexInnerBoun ) ;

MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;

MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 

MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% Verify that the internal nodes have been properly identified 
% figure
% plot( MESH.nodes(1,:) , MESH.nodes(2,:) , 'o' )
% hold on 
% plot( MESH.innerNodes(1,:) , MESH.innerNodes(2,:) , '*' ) 
% legend('Nodes' , 'Inner Nodes')

% Print some information
fprintf(' * Number of Nodes           = %d \n',MESH.numNodes);
fprintf(' * Number of inner Nodes     = %d \n',MESH.numInnerNodes);


%% 






