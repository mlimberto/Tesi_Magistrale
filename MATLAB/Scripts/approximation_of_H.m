addpath('../Inverse_Problem/')
addpath('../Inverse_Problem/Examples/')

x = -1:0.001:1 ;
h = figure ;
plot( x , 0.5 + 0.5*sign(x) , 'Linewidth' , 2 );
hold on 
tau = 0.4 ;
plot( x , smoothLS( x , tau )  , 'Linewidth' , 2 )
tau = 0.2 ;
plot( x , smoothLS( x , tau )  , 'Linewidth' , 2 )
tau = 0.1 ;
plot( x , smoothLS( x , tau )  , 'Linewidth' , 2 )
tau = 0.05 ;
plot( x , smoothLS( x , tau )  , 'Linewidth' , 2 )

legend('H(x)' , 'tau = 0.4', 'tau = 0.2', 'tau = 0.1', 'tau = 0.05')
