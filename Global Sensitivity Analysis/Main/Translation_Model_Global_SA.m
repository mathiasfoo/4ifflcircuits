%%% Global sensitivity analysis using Latin Hypercule Sampling approach.
clear all
close all
clc

LB = 10^-2;
UB = 10^2; 

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

beta_TL = 1e4;        %%% TL- Only %%%

delta_xy = 0.01;    
delta_xz = 0.01; 

%%% Number of Parameters 15 

 p = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_TL gamma_y gamma_z delta_xy delta_xz];

 ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g alpha_y_prime ...
     delta_y_prime beta_TL gamma_y gamma_z delta_xy delta_xz];

        %%%%%%% Global SA Metric Calculation %%%%%% 

N = 10000;
Point = zeros(15,N);
for ij = 1:15
    Point(ij,:) = randperm(N);
end

Interval = linspace(LB,UB,N+1);
Sample = zeros(N,15);

for i = 1:N
    for j = 1:15
        Sample(i,j) = rand*(Interval(Point(j,i)+1) - Interval(Point(j,i))) + Interval(Point(j,i));
    end
end



t_rise = zeros(N,1);
pulse_width = zeros(N,1);
pulse_height = zeros(N,1);
ss_diff = zeros(N,1);
observability = zeros(N,1);


tspan = 0:0.1:10*60^2; %%% 10 hours simulation
t_cutoff = (length(tspan)-1)/2; %%% 5 hour cutoff
ss_tolerance = 1e-5;
success_trial = []; %%% store trials that have reached steady state.
fail_trial = []; %%% store trials did not reach steady state.

Gene_X = zeros(N,2);
Gene_Y = zeros(N,2);
Gene_Z = zeros(N,2);

for niter = 1:N

    
    p(1) = ratio(1)*Sample(niter,1);
    p(2) = ratio(2)*Sample(niter,2);
    p(3) = ratio(3)*Sample(niter,3);
    
    p(4) = ratio(4)*Sample(niter,4);
    p(5) = ratio(5)*Sample(niter,5);
    p(6) = ratio(6)*Sample(niter,6);
    
    p(7) = ratio(7)*Sample(niter,7);
    p(8) = ratio(8)*Sample(niter,8);

    p(9)  = ratio(9)*Sample(niter,9);
    p(10) = ratio(10)*Sample(niter,10);
    p(11) = ratio(11)*Sample(niter,11);
    p(12)  = ratio(12)*Sample(niter,12);
    p(13) = ratio(13)*Sample(niter,13);
    p(14) = ratio(14)*Sample(niter,14);
    p(15) = ratio(15)*Sample(niter,15);
    
    

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    x0 = [0 0 0 0 0 0 0];
    [t,x] = ode23s(@TL_ODE,tspan,x0,options,p);
    L = length(tspan);
    All_x = x.*1e9; %%% convert concentration from M to nM for visualization, dont need to
    All_t = t./60;
    
    Gene_X(niter,1) = All_x(t_cutoff,1);
    Gene_Y(niter,1) = All_x(t_cutoff,2);
    Gene_Z(niter,1) = All_x(t_cutoff,3);

    Gene_X(niter,2) = All_x(end,1);
    Gene_Y(niter,2) = All_x(end,2);
    Gene_Z(niter,2) = All_x(end,3);
    
    
    label_list = ["X mRNA"; "Y mRNA"; "Z mRNA"; "X:Y Complex"; "Y Prime"; ...
        "X:Z Complex"; "GFP"]; 
    
    delta_t = diff(All_t);
    t_increment = delta_t(1); %%%% time interval in min.
    zact = All_x(:,7);
    
    %%% process the extremely small and possibly negative value, by setting
    %%% them to 0
   
    for ij = t_cutoff:length(zact)
        if zact(ij) <= 1e-15        
           zact(ij) = 0;
        end    
    end
    
    [max_zact, locmax_zact] = max(zact);
    [min_zact, locmin_zmRNA] = min(zact);
    
    if locmax_zact <= t_cutoff && zact(end) <= 0.1*max_zact

        %%% Pulse Width
        t_begin = find(zact(1:locmax_zact) <= 0.1*max_zact);
        t_end = find(zact(locmax_zact:end) <= 0.1*max_zact); 
               
         width = (t_end(1) - t_begin(end))*t_increment; %%% in unit mins
  
            %%% Rise time
            t90 = find(zact(1:locmax_zact) <= 0.9*max_zact);
            t10 = find(zact(1:locmax_zact) >= 0.1*max_zact);
            
            rise = (t90(end) - t10(1))*t_increment; 
            
            if rise <= 150 && width >= 30 %%% mins value got from biomolecules paper
                
                success_trial = [success_trial;niter]; 
                
                pulse_width(niter) = width;
                
                
                t_rise(niter) = rise;
                %%% Pulse Height
                pulse_height(niter) = max_zact - min_zact;

                %%% Steady State
                ss_diff(niter) = zact(end); %%%% New definition

                %%% Pulse Observability
                observability(niter) = pulse_height(niter)/(pulse_height(niter)+ss_diff(niter)); 
            
            else
                
                fail_trial = [fail_trial;niter];     
                pulse_width(niter) = nan; 
                t_rise(niter) = nan; 
                %%% Pulse Height
                pulse_height(niter) = nan;
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                %%% Pulse Observability
                observability(niter) = nan;
            end
            
    else
                fail_trial = [fail_trial;niter];     
                pulse_width(niter) = nan; 
                t_rise(niter) = nan; 
                %%% Pulse Height
                pulse_height(niter) = nan;
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                %%% Pulse Observability
                observability(niter) = nan;

    end
        
    
    clear All_x All_t
   
end

save('Kinetic_Translation_Global_SA_Results.mat','Sample','t_rise','pulse_width','pulse_height',...
  'observability','success_trial', 'fail_trial','ss_diff','Gene_X','Gene_Y','Gene_Z','-v7.3')


