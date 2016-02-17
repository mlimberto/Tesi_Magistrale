function h = plot_compare_results_2( c1 , c2  ) 

ax1 = gca( c1 ) ; 
ax2 = gca( c2 ) ; 

h = figure; %create new figure

temp = h.Position ;
temp(3) = 812 ;
temp(4) = 413 ;
h.Position = temp ;


s1 = subplot(1,3,1); %create and get handle to the subplot axes
s2 = subplot(1,3,2);

fig1 = get(ax1,'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
s1.XLim = [-4 4]
s1.YLim = [-4 7]
% s1.Position = [0.05 0.1 0.25 0.8]
% axis equal
s1.Position = [0.05 0.1 0.25 0.8]

fig2 = get(ax2,'children');
copyobj(fig2,s2);
s2.XLim = [-4 4]
s2.YLim = [-4 7] ;
s2.Position = [0.35 0.1 0.25 0.8]
% axis equal

colormap jet

colorbar('Position', [0.92 0.1 0.03 0.8])

% Do not display axis
axis(s1,'off')
axis(s2,'off')


end

