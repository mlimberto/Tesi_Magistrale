% Plot mesh domain
addpath('../FEM_library/')
addpath('../Inverse_Problem/')
addpath('../Inverse_Problem/Examples/')



%% Plot square

% meshFileName = '../Mesh/Square/square_coarse' ;
meshFileName = '../Mesh/Square/square' ;
% meshFileName = '../Mesh/Square/square_fine' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

h1 = figure;
pdeplot(vertices,[],elements(1:3,:));
set( h1.Children.Children , 'LineWidth' , 2 ); 
axis equal



%% Plot ovetto


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
[numElemDof,numBoundaryDof]  = select2D('P1');
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;
MESH.indexInnerElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_HEART_REGION ) ; 
MESH.indexOuterElem = find( MESH.elements( numElemDof+ 1 ,:)== DATA.FLAG_TORSO_REGION ) ; 
MESH.innerElements = MESH.elements( : , MESH.indexInnerElem ) ;
MESH.indexInnerNodes = unique( MESH.innerElements(1:numElemDof , :) ) ;
MESH.innerNodes = MESH.nodes( : , MESH.indexInnerNodes ) ; 
MESH.numInnerNodes      = size(MESH.innerNodes , 2) ;
MESH.numInnerElem       = size(MESH.innerElements,2);

%%

% Heart

h2 = figure;
ax2 = gca ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem));
axis equal

% Heart + torso

h3 = figure;
ax3 = gca ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexOuterElem));
% xlim([-7 15])
% ylim([-10 14])
axis equal
hold on 
pdeplot(MESH.vertices,[],MESH.elements(1:3,MESH.indexInnerElem));

h3.Children.Children(1).Color = [1 0 0 ];


%% Put them together in a combined plot

h4 = figure ; 
s1 = subplot(1 ,2 ,1); %create and get handle to the subplot axes
s2 = subplot(1,2,2);

fig1 = get(ax3,'children'); %get handle to all the children in the figure
% axis equal
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes


fig2 = get(ax2,'children');
% axis equal
copyobj(fig2,s2);

s1.Position = [0.05 0.05 0.5 0.9]
s1.XLim = [-7 15]

s2.Position = [0.6 0.05 0.35 0.9]










