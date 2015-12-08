% Plot mesh domain
addpath('../FEM_library/')
addpath('../Inverse_Problem/')
addpath('../Inverse_Problem/Examples/')

% meshFileName = '../Mesh/Square/square_coarse' ;
meshFileName = '../Mesh/Square/square' ;
% meshFileName = '../Mesh/Square/square_fine' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

h = figure;
pdeplot(vertices,[],elements(1:3,:));
set( h.Children.Children , 'LineWidth' , 2 ); 
axis equal
