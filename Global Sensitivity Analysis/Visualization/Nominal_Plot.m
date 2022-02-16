%%% Nominal Plot Code %%% 

clear all
close all
clc

par.P_x = 1e-9;
par.P_y = 1e-9;
par.P_z = 1e-9; 

alpha_x = 0.5;
alpha_y = 0.5;
alpha_z = 0.5;

delta_x = 0.001;
delta_y = 0.001;
delta_z = 0.001;

alpha_g = 0.5;
delta_g = 0.0005;

gamma_y = 5e3;
gamma_z = 5e3; 

alpha_y_prime = 0.5;
delta_y_prime = 0.0005;

delta_xy = 0.01;    
delta_xz = 0.01; 

beta_TX = 5e6;        %%% TX-Only %%%
beta_TL = 1e4;        %%% TL- Only %%%

beta_HY_1 = 3e4;      %%% HYBRID TY-1
beta_HY_2 = 2e6;      %%% HYBRID TY-2

nu =1e-6; 

tspan = 0:0.1:10*60^2; %%% 10 hours simulation
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

%%% TX-Only %%%%


ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...
    beta_TX gamma_y gamma_z nu];

p = ratio;


x0 = [0 0 0 0 0 0 0];
[t,x] = ode23s(@TX_ODE,tspan,x0,options,p);
All_x_TX = x(:,5).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;
clear p x t x0 Px Pz Py

%%%% TL- Only %%%% 


p = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_TL gamma_y gamma_z delta_xy delta_xz];

 
x0 = [0 0 0 0 0 0 0]; 
[t,x] = ode23s(@TL_ODE,tspan,x0,options,p);
All_x_TL = x(:,7).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;

Out = x.*10^9;
Simu_t = t./60;


clear p x t x0 Px Pz Py


%%% HY_TY1 %%% 


 p = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...  
      gamma_y gamma_z beta_HY_1 nu alpha_y_prime delta_y_prime];
  
x0 = [0 0 0 0 0 0 0];
[t,x] = ode23s(@HY_TY1_ODE_mod,tspan,x0,options,p);
All_x_HY_TY1 = x(:,4).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;

        Out = x.*10^9;
        Simu_t = t./60;


clear p x t x0 Px Pz Py 

%%% HY_TY2 %%% 

p = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_HY_2 gamma_y gamma_z delta_xy delta_xz];
 
 
x0 = [0 0 0 0 0 0 0 0];
[t,x] = ode23s(@HY_TY2_ODE,tspan,x0,options,p);
All_x_HY_TY2 = x(:,8).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;
clear p x t x0 Px Pz Py

c = 0.7; 
wid = 6; 
text = 30; 
figure
plot(All_t, All_x_TX, 'Color', [0 0.4470 0.7410 c ],'LineWidth', wid)
xlabel('Time (min)')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
xlim([0 300])
hold on
plot(All_t, All_x_TL, 'Color', [0.4660 0.6740 0.1880 c], 'LineWidth', wid)
xlabel('Time (min)')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
xlim([0 300])
hold on 
plot(All_t, All_x_HY_TY1,'Color', [0.9290, 0.6940, 0.1250, c], 'LineWidth', wid)
xlabel('Time (min)')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
xlim([0 300])
hold on
plot(All_t, All_x_HY_TY2,'Color', [0.6350 0.0780 0.1840 c], 'LineWidth', wid)
xlabel('Time (min)')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
xlim([0 300])
title('Nominal Plot')
ylabel('GFP (a.u.)')
legend('TX Circuit', 'TL Circuit', 'HY-1 Circuit', 'HY-2 Circuit') 