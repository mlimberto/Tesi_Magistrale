function h = export_plot_heart(c1)

set(c1,'Visible','on');
c1.Position(3) = 420;
c1.Position(4) = 600;
ax1 = gca(c1) ;
ax1.XLim = [-4 4]
ax1.YLim = [-4 7]
axis off
ax1.Position = [0.01 0.01 0.98 0.98]

h = c1 

end