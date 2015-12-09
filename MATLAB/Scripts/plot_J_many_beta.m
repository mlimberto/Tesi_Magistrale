cd ../../Dropbox' (Personal)'/University/PoliMi/Tesi_Magistrale/Export/

h = open('Deterministic_1/J.fig' )
J1 = h.Children.Children.YData ;
h = open('Deterministic_3/J.fig' )
J2 = h.Children.Children.YData ;
h = open('Deterministic_2/J.fig' )
J3 = h.Children.Children.YData ;
close all

h = figure
loglog( J1 , 'LineWidth' , 2 ) ; 
hold on
loglog( J2 , 'LineWidth' , 2 ) ; 
loglog( J3 , 'LineWidth' , 2 ) ; 

legend('1e-4' , '1e-6' , '1e-8') ;
grid on 

xlabel('Iteration')
ylabel('Objective function')