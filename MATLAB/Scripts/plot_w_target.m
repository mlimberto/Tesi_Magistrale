%%  Add path to libraries
addpath('../FEM_library/')
addpath('../Inverse_Problem/')
addpath('../Inverse_Problem/Examples/')

%% Plot target %% Import mesh

addpath('../Mesh/Ovetto_Circondato/') ; % add path to mesh geometry specification function

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

MESH.nodes = vertices;

% Update Mesh data with BC information and geometrical maps
[numElemDof,numBoundaryDof]  = select2D(fem);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;
MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_HEART_REGION ) ; 
MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;
MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 
MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

% Define a target control function 

tau = 0.2 ; % smoothing level
R = 1.2 ; % radius

% 1
D = [-3.9 0]; % displacement
w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;
w_target_bar_1 = extend_with_zero( w_target , MESH ) ;

% 2
D = [0 6.8]; % displacement
w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;
w_target_bar_2 = extend_with_zero( w_target , MESH ) ;

%3
D = [+3.9 0]; % displacement
w_target = circularLS( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , R , D ) ;
w_target = 1 - smoothLS(w_target , tau) ;
w_target_bar_3 = extend_with_zero( w_target , MESH ) ;

%%

c1 = figure('Visible','off') ;
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar_1(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar_1(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','off', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis([-4 4 -4 7])
    axis equal
c2 = figure('Visible','off') ;
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar_2(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar_2(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','off', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])
    axis([-4 4 -4 7])
    axis equal
c3 = figure('Visible','off') ;
    pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem),'xydata',w_target_bar_3(1:MESH.numVertices),'xystyle','interp',...
    'zdata',w_target_bar_3(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','off', 'mesh' , 'off' );
    colormap(jet);
    lighting phong
    view([0 90])    
    axis([-4 4 -4 7])
    axis equal

%% 

h1 = export_plot_heart(c1);
h2 = export_plot_heart(c2);
h3 = export_plot_heart(c3);









