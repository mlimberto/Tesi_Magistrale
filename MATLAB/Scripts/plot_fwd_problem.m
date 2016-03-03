%% Plot zd solutions

% zd shall be precomputed and stored in memory

h = figure ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
    'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off'  );
colormap(jet); lighting phong;
view([0 90])
    
% change axis
h.Children(3).XLim = [-7 15] ;
h.Children(3).YLim = [-10 14] ;

% change position in figure to reduce margins
h.Children(3).Position = [0.05 0.05 0.8 0.9] ;

% set picture position
h.Position = [446 304 554 494] ;

% equalize axes
axis equal

% slightly move colorbar

h.Children(2).Position = h.Children(2).Position + [0.03 0 0 0] ;


%% Plot difference of two solutions

% zd0 and zd1 shall be pre-computed

zd =  zd1 - zd0 ;

h = figure ;
pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',zd(1:MESH.numVertices),'xystyle','interp',...
    'zdata',zd(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on', 'mesh' , 'off'  );
colormap(jet); lighting phong;
view([0 90])
    
% change axis
h.Children(3).XLim = [-7 15] ;
h.Children(3).YLim = [-10 14] ;

% change position in figure to reduce margins
h.Children(3).Position = [0.05 0.05 0.8 0.9] ;

% set picture position
h.Position = [446 304 554 494] ;

% equalize axes
axis equal

% slightly move colorbar

h.Children(2).Position = h.Children(2).Position + [0.03 0 0 0] ;



