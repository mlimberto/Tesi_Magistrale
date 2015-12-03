function [ u ] = solveDiffusion( MESH , FE_SPACE , DATA , source_location )
    
    % Set source term according to input parameter

    if nargin<4 || isempty(source_location)
        DATA.param = [0 ; 0 ];
    else
        DATA.param = [ source_location(1) ; source_location(2) ] ;
    end
  
    % Assemble matrix and right-hand side
%     fprintf('\n Assembling ... ');
%     t_assembly = tic;
    [A, F]              =  Assembler_2D(MESH, DATA, FE_SPACE);
%     t_assembly = toc(t_assembly);
%     fprintf('done in %3.3f s', t_assembly);

    % Apply boundary conditions
%     fprintf('\n Apply boundary conditions ');
    [A_in, F_in, u_D]   =  ApplyBC_2D(A, F, FE_SPACE, MESH, DATA);


    % Solve
%     fprintf('\n Solve Au = f ... ');
%     t_solve = tic;
    u                         = zeros(MESH.numNodes,1);
    u(MESH.internal_dof)      = A_in \ F_in;
    u(MESH.Dirichlet_dof)     = u_D;
%     t_solve = toc(t_solve);
%     fprintf('done in %3.3f s \n', t_solve);

end

