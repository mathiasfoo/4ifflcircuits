%%%%%
%%%%%
%%%%% Local Sensativity Analysis
%%%%% Metrics: "Rise time, Pulse width, Pulse height, Steady state gain"

clear all
% close all
clc
global par

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

omega = 1e5; %%% TetR to Pz
gamma_y = 5e3; %%% X to Py
gamma_z = 5e3; %%% X to Pz

p = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...
    omega gamma_y gamma_z];

ratio = [alpha_x alpha_y alpha_z delta_x delta_y delta_z alpha_g delta_g ...
    omega gamma_y gamma_z];

paralist = [{'alpha_x'},{'alpha_y'},{'alpha_z'},{'delta_x'},{'delta_y'},{'delta_z'},...
    {'alpha_g'},{'delta_g'},{'omega'},{'gamma_y'},{'gamma_z'},{'nu'}];

j = 1;
idx = [];

% Nominal parameter at ii = 28

for k = 1:12 %%% Including a Loop for each Variable %%%
    
    for ii = [0.001:0.001:0.009 0.01:0.01:0.09 0.1:0.1:0.9 1:1:9 10:10:90 100:100:1000]
        j
        if k == 1   %%% Initialize the Previous Variable
            p(k) = ii*ratio(k);
        else
            p(k-1) = ratio(k-1);
            p(k) = ii*ratio(k);
        end
        
        tspan = 0:0.1:10*60^2;          %%% 10 hours simulation
        t_cutoff = (length(tspan)-1)/2; %%% 5 hour cutoff
        ss_tolerance = 1e-5;
        ss_trial = [];                  %%% store trials that have reached steady state.
        nss_trial = [];                 %%% store trials that have not reached steady state.
        
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        x0 = [0 0 0 0 0];
        
        [t,x] = ode23s(@TX_ODE,tspan,x0,options,p);
        L = length(tspan);
        All_x = x.*1e0;                  %%% convert from M to nM for visualization
        All_t = t./60;
        
        figure(1)
        xlabel('Time (min)')
        ylabel('GFP Output')
        plot(All_t,All_x(:,5)*1e9,'LineWidth',2)
        hold on
        title(['\',paralist{k}])
        
        %%% Calculating metric for GFP
        
        delta_t = diff(All_t);
        gfp = All_x(:,5)./1e0;
        de_avg = abs(gfp(t_cutoff) - gfp(end));
        
        [max_gfp, locmax_gfp] = max(gfp);
        [min_gfp, locmin_gfp] = min(gfp);
        
        if locmax_gfp*delta_t(1) <= t_cutoff*delta_t(1) && gfp(end) <= 0.1*max_gfp
            
            [max_gfp, locmax_gfp] = max(gfp);
            [min_gfp, locmin_gfp] = min(gfp);
            t90 = find(gfp(1:locmax_gfp) <= 0.9*max_gfp);
            t10 = find(gfp(1:locmax_gfp) >= 0.1*max_gfp);
            
            
            %%% Calculating Metrics: "Rise time, Pulse width, Pulse height, Steady state gain"
            
            %%% Rise time
            t_rise(k,j) = (t90(end) - t10(1))*delta_t(1);
            
            %%% Pulse Width          
            t_begin = find(gfp(1:locmax_gfp) <= 0.1*max_gfp);
            t_end = find(gfp(locmax_gfp:end) <= 0.1*max_gfp);
            if isempty(t_end) == 1
                pulse_width(k,j) = nan;
            else
                pulse_width(k,j) = (t_end(1) - t_begin(end))*delta_t(1);
            end
            
            %%% Pulse Height
            pulse_height(k,j) = max_gfp - min_gfp;
            
            %%% Steady State Gain
            ss_gain(k,j) = gfp(end);  %%% Inital Concentration is equal to zero
        else
            nss_trial = [nss_trial;j];
            
            %%% Rise Time
            t_rise(k,j) = nan;
            %%% Pulse Width
            pulse_width(k,j) = nan;
            %%% Pulse Height
            pulse_height(k,j) = nan;
            %%% Steady State Gain
            ss_gain(k,j) = nan;
            
        end
        
        idx = [idx ii];
        All_x_idx(j,:,:) = All_x(:,5);
        [t_rise(k,j) pulse_width(k,j) pulse_height(k,j) ss_gain(k,j)];
        
        j = j + 1;
        if j > 55
            j = 1;
        end
    end
end

% save txmetricnoscale t_rise pulse_width pulse_height ss_gain

%% Plotting the Metrics
% figure(2)
% subplot(2,2,1)
% semilogx(idx,t_rise,'o-','LineWidth',2)
% title('Rise Time')
% set(gca,'FontName', 'Times New Roman')
% ylim([0 0.3])
% xlim([1e-3 1e3])
%
% subplot(2,2,2)
% semilogx(idx,pulse_width,'o-','LineWidth',2)
% title('Pulse Width')
% set(gca,'FontName', 'Times New Roman')
% ylim([0 0.3])
% xlim([1e-3 1e3])
%
% subplot(2,2,3)
% semilogx(idx,pulse_height,'o-','LineWidth',2)
% title('Pulse Height')
% set(gca,'FontName', 'Timew New Roman')
% ylim([0 0.3])
% xlim([1e-3 1e3])
%
% subplot(2,2,4)
% semilogx(idx,ss_gain,'o-','LineWidth',2)
% title('Seady State Gain')
% set(gca,'FontName', 'Timew New Roman')
% ylim([0 0.3])
% xlim([1e-3 1e3])

%%% Need to Save Results %%%
