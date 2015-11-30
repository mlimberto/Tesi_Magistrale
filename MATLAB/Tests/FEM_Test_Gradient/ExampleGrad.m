%% Import adjoint problem example

cd ../FEM_Test_Adj_Problem/

ExampleAdj ;

cd ../FEM_Test_Gradient/

% We now want to evaluate the Riesz element corresponding to the B'*p
% operator which maps the control function to the real numbers

%% Assemble lhs matrices 

% First we assemble the reaction matrix in the heart subdomain

fprintf('\n Assembling reaction matrix in the heart subdomain ... ');
t_assembly_react = tic;
A_react = Assembler_2D( MESH , DATA , FE_SPACE , 'reaction' , [] , [] , FLAG_HEART_REGION ) ; 
t_assembly_react = toc(t_assembly_react);
fprintf('done in %3.3f s\n', t_assembly_react);
clear t_assembly_react ;

% We now assemble the Laplacian matrix in the heart subdomain 

% REMARK this matrix has already been assembled and goes under the name
% A_source


% Eventually we define the lhs gradient matrix 

A_grad = A_react + A_source ; 


%% Assemble rhs matrix

% We wish to assemble the matrix evaluating the source term B' * p 

% We can obtain it again from the A_source matrix
A_tBp = - DATA.coeffRhs * A_source ; 

%% Evaluate rhs 

F_grad = A_tBp * p ; 

% Visualize rhs 
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',F_grad(1:MESH.numVertices),'xystyle','interp',...
       'zdata',F_grad(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off' );
colormap(jet);
lighting phong
title('Gradient source term')

%% Restrict only to the heart subdomain 

A_grad_restricted = A_grad( MESH.indexInnerNodes , MESH.indexInnerNodes ) ; 

F_grad_restricted = F_grad( MESH.indexInnerNodes ) ; 

%% Solve the problem 

dw = A_grad_restricted \ F_grad_restricted ;  

%% Visualize dw

dwbar = extend_with_zero( dw , MESH ) ; 
figure
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',dwbar(1:MESH.numVertices),'xystyle','interp',...
       'zdata',dwbar(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on', 'mesh' , 'off' );
colormap(jet);
lighting phong
title('Control derivative Riesz element')


%% Update w 

% Apply gradient step
w_new = ( 1 - DATA.beta*DATA.gstep )*w - DATA.gstep * dw ; 

% Visualize it
H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w_new ) ; 
[X,Y] = meshgrid(-1:0.02:1) ; 
figure
surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp ; colormap jet ; title('Control function') ; axis equal ;
clear H ; clear X ; clear Y ; 

% Apply projection step 
w_new = min( max( w_new , 0 ) , 1 )  ;

% Visualize it
H = scatteredInterpolant( MESH.innerNodes(1,:)' , MESH.innerNodes(2,:)' , w_new ) ; 
[X,Y] = meshgrid(-1:0.02:1) ; 
figure
surf(X,Y,H(X,Y) , 'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp ; colormap jet ; title('Updated control function') ; axis equal ;
clear H ; clear X ; clear Y ; 
