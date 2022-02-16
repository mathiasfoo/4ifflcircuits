%%%% Global Sensativity Analysis Results
%%%% Violin Plot Code 
close all
clear all
clc

%%% Global Results from Global-SA TX Circuit %%%

load Kinetic_Transcription_Global_SA_Results.mat;
RT_TX = t_rise; 
PH_TX = pulse_height;
PW_TX = pulse_width; 
SSG_TX = ss_diff;
PQ_TX = observability;

PE_TX = pulse_width.*pulse_height; 

A_TX = (pulse_width)./(t_rise+pulse_width);

Passed_TX = length(success_trial);
Faild_TX = length(fail_trial);
clear t_rise pulse_height pulse_width ss_diff fail_trial observability success_trial

%%% Global Results from Global-SA TL Circuit %%%

load Kinetic_Translation_Global_SA_Results.mat;
RT_TL = t_rise; 
PH_TL = pulse_height;
PW_TL = pulse_width; 
SSG_TL = ss_diff;
PQ_TL = observability;

A_TL = (pulse_width)./(t_rise+pulse_width);

PE_TL = pulse_width.*pulse_height; 

Passed_TL = length(success_trial);
Faild_TL = length(fail_trial);
clear t_rise pulse_height pulse_width ss_diff fail_trial observability success_trial


%%% Loal Results from Global-SA HY-TY1 Modified Circuit %%%

load Kinetic_Hybrid_TY1_mod_Global_SA_Results.mat;
RT_HY_TY1_mod = t_rise; 
PH_HY_TY1_mod = pulse_height;
PW_HY_TY1_mod = pulse_width; 
SSG_HY_TY1_mod = ss_diff;
PQ_HY_TY1_mod = observability;

A_HY_TY1_mod = (pulse_width)./(t_rise+pulse_width);

PE_HY_TY1_mod = pulse_width.*pulse_height; 

Passed_HY_TY1_mod = length(success_trial);
Faild_HY_TY1_mod = length(fail_trial);
clear t_rise pulse_height pulse_width ss_diff fail_trial observability success_trial

%%% Loal Results from Global-SA HY_TY2 Circuit %%%

load Kinetic_Hybrid_TY2_Global_SA_Results.mat;
RT_HY_TY2 = t_rise; 
PH_HY_TY2 = pulse_height;
PW_HY_TY2 = pulse_width; 
SSG_HY_TY2 = ss_diff;
PQ_HY_TY2 = observability;

PE_HY_TY2 = pulse_width.*pulse_height; 

A_HY_TY2 = (pulse_width)./(t_rise+pulse_width);

Passed_HY_TY2 = length(success_trial);
Faild_HY_TY2 = length(fail_trial);
clear t_rise pulse_height pulse_width ss_diff fail_trial observability success_trial


%%%%%%%%%%%%%%% Visulization %%%%%%%%%%%%%%%%%%

%%% Swarm Scatter Chart
N = 10000;
names = ["TX","TL","Hybrid TY-1", "Hybrid TY-2"];
x = categorical([ones(1,N) 2*ones(1,N) 3*ones(1,N) 4*ones(1,N)],[1 2 3 4], names);

name = {'TX', 'TL', 'Hybrid TY-1', 'Hybrid TY-2'};

text = 20; 

figure
subplot(2,2,1) 
vs = Violin((RT_TX), 1);
vs = Violin((RT_TL), 2);
vs = Violin((RT_HY_TY1_mod),3); 
vs = Violin((RT_HY_TY2),4); 
ylabel('Rise Time [min]');
title('Rise Time')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})

subplot(2,2,2) 
vs = Violin(PW_TX, 1);
vs = Violin(PW_TL, 2);
vs = Violin(PW_HY_TY1_mod,3); 
vs = Violin(PW_HY_TY2,4); 
ylabel('Pulse Width [min]');
title('Pulse Width')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text)
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})

subplot(2,2,3) 

vs = Violin(log(PH_TX), 1, 'Width', 0.45);
vs = Violin(log(PH_TL), 2);
vs = Violin(log(PH_HY_TY1_mod),3); 
vs = Violin(log(PH_HY_TY2),4, 'Width', 0.45); 
ylabel('log(Pulse Height)');
title('Pulse Height')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})

%%%%%%Modification
SSG_TX_logabs = log(abs(SSG_TX));
SSG_TL_logabs = log(abs(SSG_TL));
SSG_TY1_logabs = log(abs(SSG_HY_TY1_mod));
SSG_TY2_logabs = log(abs(SSG_HY_TY2));
txj= 0;
tlj= 0;
ty1j= 0;
ty2j= 0;

for i = 1:10000
    if ~isnan(SSG_TX_logabs(i)) && SSG_TX_logabs(i)>=-10^100  %%%% Process the data to get rid of -inf and NaN
        txj=txj+1;
        CTX(txj)=SSG_TX_logabs(i);
    end
    
    if ~isnan(SSG_TL_logabs(i)) && SSG_TL_logabs(i)>=-10^100  %%%% Process the data to get rid of -inf and NaN
        tlj=tlj+1;
        CTL(tlj)=SSG_TL_logabs(i);
    end
    
    if ~isnan(SSG_TY1_logabs(i)) && SSG_TY1_logabs(i)>=-10^100  %%%% Process the data to get rid of -inf and NaN
        ty1j=ty1j+1;
        CTY1(ty1j)=SSG_TY1_logabs(i);
    end
    
    if ~isnan(SSG_TY2_logabs(i)) && SSG_TY2_logabs(i)>=-10^100  %%%% Process the data to get rid of -inf and NaN
        ty2j=ty2j+1;
        CTY2(ty2j)=SSG_TY2_logabs(i);
    end
    
end


subplot(2,2,4) 
vs = Violin(CTX, 1);
vs = Violin(CTL, 2);
vs = Violin(CTY1,3); 
vs = Violin(CTY2,4); 
ylabel('log(|Steady State Value|)');
title('Steady State')
set(gca,'FontName', 'Times')
set(gca,'FontSize', text) 
set(gca,'XTick',[])
xticks([1 2 3 4])
xticklabels({'TX', 'TL', 'HY-1', 'HY-2'})
legend('Simulation Data', 'Kernal Density Estimate', 'Interquartile Range (IQR)', '', '','Median'); 


