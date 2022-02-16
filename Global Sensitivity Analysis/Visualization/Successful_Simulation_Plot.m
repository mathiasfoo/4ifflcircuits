%%%% Plot Simulation %%%%

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

beta_HY_1 = 1e5;     %%% HYBRID TY-1
beta_HY_2 = 2e6;     %%% HYBRID TY-2

nu =1e-6; 

num = 25; %%% number of samples

%c = [0.1 0.25 0];
tspan = 0:0.1:10*60^2; %%% 10 hours simulation
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
figure 

%%% Global Results from Global-SA TX Circuit %%%

load Kinetic_Transcription_Global_SA_Results.mat;
Passed_TX = length(success_trial);
TX_Selection = randi(Passed_TX, 1, num); 

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...
    beta_TX gamma_y gamma_z nu];

for i = 1:num
for j = 1:12    
    
p(j) = ratio(j)*Sample(success_trial(TX_Selection(i)),j);

end

x0 = [0 0 0 0 0 0 0];
[t,x] = ode23s(@TX_ODE,tspan,x0,options,p);
All_x_TX = x(:,5).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;


hold on
subplot(2,2,1)
plot(All_t, All_x_TX, 'LineWidth', 2)
title ('TX Circuit') 
xlabel('Time (min)')
xlim([-10 300])
set(gca,'FontName', 'Times')
set(gca,'FontSize', 18)


end
hold off
clear success_trial Sample p

%% Global Results from Global-SA TL Circuit %%%
load Kinetic_Translation_Global_SA_Results.mat;
Passed_TL = length(success_trial);
TL_Selection = randi(Passed_TL, 1, num);

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_TL gamma_y gamma_z delta_xy delta_xz];

for i = 1:num
for j = 1:15  
    
p(j) = ratio(j)*Sample(success_trial(TL_Selection(i)),j);

end 

 
x0 = [0 0 0 0 0 0 0]; 
[t,x] = ode23s(@TL_ODE,tspan,x0,options,p);
All_x_TL = x(:,7).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;

hold on
subplot(2,2,2)
plot(All_t, All_x_TL, 'LineWidth', 2)
title ('TL Circuit') 
xlabel('Time (min)')
xlim([-10 300])
set(gca,'FontName', 'Times')
set(gca,'FontSize', 18)


end
hold off
clear success_trial Sample p


  
%% Loal Results from Global-SA HY-TY1 Modified Circuit %%%
load Kinetic_Hybrid_TY1_mod_Global_SA_Results.mat;
Passed_HY_TY1_mod = length(success_trial);
HY_TY1_mod_Selection = randi(Passed_HY_TY1_mod, 1, num);

beta_HY_1 = 3e4;    

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...  
      gamma_y gamma_z beta_HY_1 nu alpha_y_prime delta_y_prime];
  
for i = 1:num
for j = 1:14 
    
p(j) = ratio(j)*Sample(success_trial(HY_TY1_mod_Selection(i)),j);

end 
 
  
x0 = [0 0 0 0 0 0 0];
[t,x] = ode23s(@HY_TY1_ODE_mod,tspan,x0,options,p);
All_x_HY_TY1_mod = x(:,4).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;


hold on
subplot(2,2,3)
plot(All_t, All_x_HY_TY1_mod, 'LineWidth', 2)
title ('HY-1 Circuit') 
xlabel('Time (min)')
xlim([-10 300])
set(gca,'FontName', 'Times')
set(gca,'FontSize', 18)

end
hold off
clear success_trial Sample p

%% Loal Results from Global-SA HY_TY2 Circuit %%%
load Kinetic_Hybrid_TY2_Global_SA_Results.mat;
Passed_HY_TY2 = length(success_trial);
HY_TY2_Selection = randi(Passed_HY_TY2, 1, num);

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_HY_2 gamma_y gamma_z delta_xy delta_xz];

for i = 1:num
for j = 1:15
    
p(j) = ratio(j)*Sample(success_trial(HY_TY2_Selection(i)),j);

end 
 
x0 = [0 0 0 0 0 0 0 0];
[t,x] = ode23s(@HY_TY2_ODE,tspan,x0,options,p);
All_x_HY_TY2 = x(:,8).*1e9; %%% convert concentration from M to nM for visualization, dont need to
All_t = t./60;

hold on
subplot(2,2,4)
plot(All_t, All_x_HY_TY2,  'LineWidth', 2)
title ('HY-2 Circuit') 
xlabel('Time (min)')
xlim([-10 300])
set(gca,'FontName', 'Times')
set(gca,'FontSize', 18)

end
hold off
clear success_trial Sample p

