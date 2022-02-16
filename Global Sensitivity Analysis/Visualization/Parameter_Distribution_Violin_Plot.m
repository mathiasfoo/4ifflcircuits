%%%% Plot Simulation %%%%
close all
clear all
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

beta_HY_1 = 3e4;     %%% HYBRID TY-1
beta_HY_2 = 2e6;     %%% HYBRID TY-2

nu =1e-6; 


%%% Global Results from Global-SA TX Circuit %%%

load Kinetic_Transcription_Global_SA_Results.mat;
Passed_TX = length(success_trial);

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...
    beta_TX gamma_y gamma_z nu];

for i = 1:Passed_TX
for j = 1:12    
    
p(i,j) = ratio(j)*Sample(success_trial(i),j);

end
    alpha_x_TX(i) = p(i,1);  
    alpha_y_TX(i) = p(i,2); 
    alpha_z_TX(i) = p(i,3); 
    delta_x_TX(i) = p(i,4);
    delta_y_TX(i) = p(i,5); 
    delta_z_TX(i) = p(i,6); 
    alpha_g_TX(i) = p(i,7); 
    delta_g_TX(i) = p(i,8);   
    beta_TX(i) = (p(i,9)-ratio(9))/ratio(9);
    gamma_y_TX(i) = p(i,10); 
    gamma_z_TX(i) = p(i,11); 
    nu_TX(i) = p(i,12); 
end

par_TX = p;

clear success_trial Sample p

%%% Global Results from Global-SA TL Circuit %%%
load Kinetic_Translation_Global_SA_Results.mat;
Passed_TL = length(success_trial);
TL_Selection = randi(Passed_TL, 1, 10);

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_TL gamma_y gamma_z delta_xy delta_xz];

for i = 1:Passed_TL
for j = 1:15    
p(i,j) = ratio(j)*Sample(success_trial(i),j);

end
    alpha_x_TL(i) = p(i,1);  
    alpha_y_TL(i) = p(i,2); 
    alpha_z_TL(i) = p(i,3); 
    delta_x_TL(i) = p(i,4);
    delta_y_TL(i) = p(i,5); 
    delta_z_TL(i) = p(i,6); 
    alpha_g_TL(i) = p(i,7); 
    delta_g_TL(i) = p(i,8); 
    alpha_y_prime_TL(i) = p(i,9); 
    delta_y_prime_TL(i) = p(i,10); 
    beta_TL (i) = (p(i,11)-ratio(11))/ratio(11); 
    gamma_y_TL (i) = p(i,12); 
    gamma_z_TL (i) = p(i,13); 
    delta_xy_TL (i) = p(i,14); 
    delta_xz_TL(i) = p(i,15); 

end 

par_TL = p;

clear success_trial Sample p

%%% Loal Results from Global-SA HY-TY1 Modified Circuit %%%
load Kinetic_Hybrid_TY1_mod_Global_SA_Results.mat;
Passed_HY_TY1_mod = length(success_trial);
   

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...  
      gamma_y gamma_z beta_HY_1 nu alpha_y_prime delta_y_prime];
  
for i = 1:Passed_HY_TY1_mod
for j = 1:14 
    
p(i, j) = ratio(j)*Sample(success_trial(i),j);

end 
    alpha_x_HY_TY1_mod(i) = p(i,1);  
    alpha_y_HY_TY1_mod(i) = p(i,2); 
    alpha_z_HY_TY1_mod(i) = p(i,3); 
    delta_x_HY_TY1_mod(i) = p(i,4);
    delta_y_HY_TY1_mod(i) = p(i,5); 
    delta_z_HY_TY1_mod(i) = p(i,6); 
    alpha_g_HY_TY1_mod(i) = p(i,7); 
    delta_g_HY_TY1_mod(i) = p(i,8);   
    gamma_y_HY_TY1_mod(i) = p(i,9); 
    gamma_z_HY_TY1_mod(i) = p(i,10); 
    beta_HY_TY1_mod (i) = (p(i,11)-ratio(11))/ratio(11);
    nu_HY_TY1_mod (i) = p(i,12); 
    alpha_y_prime_HY_TY1_mod (i) = p(i,13); 
    delta_y_prime_HY_TY1_mod (i) = p(i,14); 

end

par_HY_TY1_mod = p; 

clear success_trial Sample p

%%% Loal Results from Global-SA HY_TY2 Circuit %%%
load Kinetic_Hybrid_TY2_Global_SA_Results.mat;
Passed_HY_TY2 = length(success_trial);

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_HY_2 gamma_y gamma_z delta_xy delta_xz];

for i = 1:Passed_HY_TY2
    for j = 1:15
    
    p(i,j) = (ratio(j)*Sample(success_trial(i),j));
 
    end 
    
    alpha_x_HY_TY2(i) = p(i,1);  
    alpha_y_HY_TY2(i) = p(i,2); 
    alpha_z_HY_TY2(i) = p(i,3); 
    delta_x_HY_TY2(i) = p(i,4);
    delta_y_HY_TY2(i) = p(i,5); 
    delta_z_HY_TY2(i) = p(i,6); 
    alpha_g_HY_TY2(i) = p(i,7); 
    delta_g_HY_TY2(i) = p(i,8); 
    alpha_y_prime_HY_TY2(i) = p(i,9); 
    delta_y_prime_HY_TY2(i) = p(i,10); 
    beta_HY_TY2 (i) = (p(i,11)-ratio(11))/ratio(11); 
    gamma_y_HY_TY2 (i) = p(i,12); 
    gamma_z_HY_TY2 (i) = p(i,13); 
    delta_xy_HY_TY2 (i) = p(i,14); 
    delta_xz_HY_TY2(i) = p(i,15); 
end

labels = {'\alpha_X','\alpha_Y','\alpha_Z','\delta_X','\delta_Y','\delta_Z',...
    '\alpha_G','\delta_G','\alpha_{Y''}','\delta_{Y''}',...
    '\beta_{HY2}','\gamma_Y','\gamma_Z','\delta_{XY}','\delta_{XZ}'};

paralist = [{'alpha_X'},{'alpha_Y'},{'alpha_Z'},{'delta_X'},{'delta_Y'},{'delta_Z'},...
    {'alpha_G'},{'delta_G'},{'alpha_{Y''}'},{'delta_{Y''}'},...
    {'beta_{HY2}'},{'gamma_Y'},{'gamma_Z'},{'delta_{XY}'},{'delta_{XZ}'}];

par_HY_TY2 = p; 
clear success_trial Sample p

%%%% Order of Dis TX,TL,HY-TY1,HY-TY1mod, HY-TY2 %%%%

%%% Swarm Scatter Chart
names = ["TX-only","TL-only","Hybrid TY-1", "Hybrid TY-2"];
x = categorical([ones(1,Passed_TX) 2*ones(1,Passed_TL) 3*ones(1,Passed_HY_TY1_mod) 4*ones(1,Passed_HY_TY2)],[1 2 3 4], names);
text = 12; 

figure 
subplot (4,3,1)    
vs = Violin((alpha_x_TX), 1);
vs = Violin((alpha_x_TL), 2);
vs = Violin((alpha_x_HY_TY1_mod),3); 
vs = Violin((alpha_x_HY_TY2),4); 
ylim([0 75])
title('\alpha_X Distribution')
ylabel('\alpha_X Value') 
set(gca,'FontSize', text) 
set(gca,'FontName', 'Times')
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on 

subplot (4,3,2)    
vs = Violin((alpha_y_TX), 1);
vs = Violin((alpha_y_TL), 2);
vs = Violin((alpha_y_HY_TY1_mod),3); 
vs = Violin((alpha_y_HY_TY2),4); 
title('\alpha_Y Distribution')
ylabel('\alpha_Y Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on 
% 
subplot (4,3,3)    
vs = Violin((alpha_z_TX), 1);
vs = Violin((alpha_z_TL), 2);
vs = Violin((alpha_z_HY_TY1_mod),3); 
vs = Violin((alpha_z_HY_TY2),4); 
title('\alpha_Z Distribution')
ylabel('\alpha_Z Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on 
% 
subplot (4,3,4)    
vs = Violin((delta_x_TX), 1);
vs = Violin((delta_x_TL), 2);
vs = Violin((delta_x_HY_TY1_mod),3); 
vs = Violin((delta_x_HY_TY2),4); 
title('\delta_X Distribution')
ylabel('\delta_X Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on 
% 
subplot (4,3,5)    
vs = Violin((delta_y_TX), 1);
vs = Violin((delta_y_TL), 2);
vs = Violin((delta_y_HY_TY1_mod),3); 
vs = Violin((delta_y_HY_TY2),4); 

title('\delta_Y Distribution')
ylabel('\delta_Y Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on   
   
subplot (4,3,6)    
vs = Violin((delta_z_TX), 1);
vs = Violin((delta_z_TL), 2);
vs = Violin((delta_z_HY_TY1_mod),3); 
vs = Violin((delta_z_HY_TY2),4); 

title('\delta_Z Distribution')
ylabel('\delta_Z Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on  

subplot (4,3,7)    
vs = Violin((alpha_g_TX), 1);
vs = Violin((alpha_g_TL), 2);
vs = Violin((alpha_g_HY_TY1_mod),3); 
vs = Violin((alpha_g_HY_TY2),4); 
title('\alpha_G Distribution')
ylabel('\alpha_G Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on   

subplot (4,3,8)    
vs = Violin((delta_g_TX), 1);
vs = Violin((delta_g_TL), 2);
vs = Violin((delta_g_HY_TY1_mod),3); 
vs = Violin((delta_g_HY_TY2),4); 

title('\delta_G Distribution')
ylabel('\delta_G Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on 

subplot (4,3,10)    
vs = Violin((gamma_y_TX), 1);
vs = Violin((gamma_y_TL), 2);
vs = Violin((gamma_y_HY_TY1_mod),3); 
vs = Violin((gamma_y_HY_TY2),4); 

title('\gamma_Y Distribution')
ylabel('\gamma_Y Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on  

subplot (4,3,11)    
vs = Violin((gamma_z_TX), 1);
vs = Violin((gamma_z_TL), 2);
vs = Violin((gamma_z_HY_TY1_mod),3); 
vs = Violin((gamma_z_HY_TY2),4);
colormap(gca, 'jet')
title('\gamma_Z Distribution')
ylabel('\gamma_Z Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on  

%%% Figure 2 %%% 


names_1 = ["TL-only", "HY-2"];
x_1 = categorical([1*ones(1,Passed_TL) 2*ones(1,Passed_HY_TY2)],[1 2], names_1);

names_2 = ["TX-only","HY-1"];
x_2 = categorical([ones(1,Passed_TX) 2*ones(1,Passed_HY_TY1_mod)],[1 2], names_2);

names_3 = ["TL-only","HY-1", "HY-2"];
x_3 = categorical([ones(1,Passed_TL) 2*ones(1,Passed_HY_TY1_mod) 3*ones(1,Passed_HY_TY2)],[1 2 3], names_3);

figure 
subplot (2,3,1)    
vs = Violin((delta_xy_TL), 1);
vs = Violin((delta_xy_HY_TY2), 2);
title('\delta X:Y Distribution')
ylabel('\delta X:Y Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2])
xticklabels({'TL','HY-2'})
hold on 

subplot (2,3,2)    
vs = Violin((delta_xz_TL), 1);
vs = Violin((delta_xz_HY_TY2), 2);
title('\delta X:Z Distribution')
ylabel('\delta X:Z Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2])
xticklabels({'TL','HY-2'})
hold on 

subplot (2,3,3)    
vs = Violin((nu_TX), 1);
vs = Violin((nu_HY_TY1_mod), 2);
title('\nu Distribution')
ylabel('\nu Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2])
xticklabels({'TX','HY-1'})
hold on 

subplot (2,3,4)    
s15.data = [alpha_y_prime_TL, alpha_y_prime_HY_TY1_mod, alpha_y_prime_HY_TY2];
vs = Violin((alpha_y_prime_TL), 1);
vs = Violin((alpha_y_prime_HY_TY1_mod), 2);
vs = Violin((alpha_y_prime_HY_TY2),3); 
title('\alpha_yprime Distribution')
ylabel('\alpha_yprime Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3])
xticklabels({'TL','HY-1', 'HY-2'})
hold on 

subplot (2,3,5)    
vs = Violin((delta_y_prime_TL), 1);
vs = Violin((delta_y_prime_HY_TY1_mod), 2);
vs = Violin((delta_y_prime_HY_TY2),3); 
title('\alpha_yprime Distribution')
title('\delta_yprime Distribution')
ylabel('\delta_yprime Value') 
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3])
xticklabels({'TL','HY-1', 'HY-2'})
hold on 

subplot (2,3,6)    
s17.data = [beta_TX, beta_TL,beta_HY_TY1_mod, beta_HY_TY2];
vs = Violin((beta_TX), 1);
vs = Violin((beta_TL),2);
vs = Violin((beta_HY_TY1_mod),3); 
vs = Violin((beta_HY_TY2),4);
title('\beta Distribution')
ylabel('\beta(New Value - Old Value)/Old Value') 
set(gca,'FontName', 'Times') 
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
hold on  
   
   
   
