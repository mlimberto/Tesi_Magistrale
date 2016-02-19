%% plot results

addpath('~/Developer/Tesi_Magistrale/MATLAB/Scripts/')
cd ~/Dropbox/University/PoliMi/Tesi_Magistrale/Export/


%% linear vs lin-nonlin vs nonline for test case 1

cd Ovetto_Deterministic/Linear_1

c1 = openfig('w_opt.fig','reuse')
% set(c1,'Visible','off')
h1 = export_plot_heart( c1 );

cd ../Linear_NonLinear_zd_1/
c2 = openfig('w_opt_fixed.fig','reuse')
set(c2,'Visible','off')
h2 = export_plot_heart( c2 );

cd ../Non_Linear_1/
c3 = openfig('w_opt.fig','reuse')
set(c3,'Visible','off')
h3 = export_plot_heart( c3 );

%% linear vs lin-nonlin vs nonline for test case 2

cd ..//Linear_2

c1 = openfig('w_opt.fig','reuse')
% set(c1,'Visible','off')
h1 = export_plot_heart( c1 );

cd ../Linear_NonLinear_zd_2/
c2 = openfig('w_opt.fig','reuse')
set(c2,'Visible','off')
h2 = export_plot_heart( c2 );

cd ../Non_Linear_2/
c3 = openfig('w_opt.fig','reuse')
set(c3,'Visible','off')
h3 = export_plot_heart( c3 );


%% linear vs lin-nonlin vs nonline for test case 3

cd ..//Linear_3

c1 = openfig('w_opt.fig','reuse')
% set(c1,'Visible','off')
h1 = export_plot_heart( c1 );

cd ../Linear_NonLinear_zd_3/
c2 = openfig('w_opt_fixed.fig','reuse')
set(c2,'Visible','off')
h2 = export_plot_heart( c2 );

cd ../Non_Linear_3/
c3 = openfig('w_opt.fig','reuse')
set(c3,'Visible','off')
h3 = export_plot_heart( c3 );



