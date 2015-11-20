%EXAMPLE2D shows how to solve a problem with the FEM library.

clc
% clear all

% Add path to FE solver
addpath('../../FEM_library/')

% Define which elements to use
fem = 'P1';

%% Import mesh

meshFileName = '../../Mesh/Square/square_coarse' ;
%meshFileName = '../../Mesh/Square/square' ;

[vertices, boundaries, elements] = msh_to_Mmesh(meshFileName, 2) ;

%% Visualize mesh

pdeplot(vertices,[],elements(1:3,:))
axis equal

%% Assemble matrix and solve problem

[U, FE_SPACE, MESH, DATA, errorL2, errorH1] = Elliptic2D_Solver(elements, vertices, boundaries, fem, 'DirichletNeumann_data');

%% Visualize solution 

figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',U(1:MESH.numVertices),'xystyle','interp',...
       'zdata',U(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on','mesh','on');
%colormap(jet);
%lighting phong





