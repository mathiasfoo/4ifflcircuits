clc
clear all
% close all


labels = {'\alpha_X','\alpha_Y','\alpha_Z','\delta_X','\delta_Y','\delta_Z',...
    '\alpha_G','\delta_G','\gamma_Y','\gamma_Z','\beta_{HY1}','\nu','\alpha_{Y''}','\delta_{Y''}'};

paralist = [{'alpha_X'},{'alpha_Y'},{'alpha_Z'},{'delta_X'},{'delta_Y'},{'delta_Z'},...
    {'alpha_G'},{'delta_G'},{'gamma_Y'},{'gamma_Z'},{'gamma_Y'},{'gamma_Z'},{'beta_{HY1}'},{'nu'},{'alpha_{Y''}'},{'delta_{Y''}'}];


% Plotting Phase Plot
load hy1metricupdate.mat

t_rise(t_rise >= 150) = nan;
m1 = isnan(t_rise);

pulse_width(pulse_width <= 30) = nan;
m2 = isnan(pulse_width);

pulse_height = pulse_height.*1e3; % -> convert to pM
pulse_height = log10(pulse_height);
m3 = isnan(pulse_height);

ss_gain(ss_gain./1e9 <= 1e-15) = nan;
ss_gain = ss_gain.*1e6; % -> convert to fM
ss_gain = log10(ss_gain);
m4 = isnan(ss_gain);

m_all = m1 | m2 | m3 | m4; 
m_all = ~m_all;

t_rise = t_rise.*m_all;
t_rise(t_rise == 0) = nan;
pulse_width = pulse_width.*m_all;
pulse_width(pulse_width == 0) = nan;
pulse_height = pulse_height.*m_all;
pulse_height(pulse_height == 0) = nan;
ss_gain = ss_gain.*m_all;
ss_gain(ss_gain == 0) = nan;

for k = 1:14
    %%% Rise time
    mintrise(k) = min(t_rise(k,:));
    maxtrise(k) = max(t_rise(k,:));
    nominaltrise(k) = t_rise(k,28);
    
    %%% Pulse width
    minpulsewidth(k) = min(pulse_width(k,:));
    maxpulsewidth(k) = max(pulse_width(k,:));
    nominalpulsewidth(k) = pulse_width(k,28);
    
    %%% Pulse height
    minpulseheight(k) = min(pulse_height(k,:)).*1e0;
    maxpulseheight(k) = max(pulse_height(k,:)).*1e0;
    nominalpulseheight(k) = pulse_height(k,28).*1e0;
    
    %%% Steady State Gain
    minssgain(k) = min(ss_gain(k,:))./1e-0;
    maxssgain(k) = max(ss_gain(k,:))./1e-0;
    nominalssgain(k) = ss_gain(k,28)./1e-0;
    
end

% Determine the maximum and minimum value and the respective parameters
[trmax, trlmax] = sort(maxtrise); disp([trmax(end), labels(trlmax(end))])
[trmin, trlmin] = sort(mintrise); disp([trmin(1), labels(trlmin(1))])

[pwmax, pwlmax] = sort(maxpulsewidth); disp([pwmax(end), labels(pwlmax(end))])
[pwmin, pwlmin] = sort(minpulsewidth); disp([pwmin(1), labels(pwlmin(1))])

[phmax, phlmax] = sort(maxpulseheight); disp([phmax(end), labels(phlmax(end))])
[phmin, phlmin] = sort(minpulseheight); disp([phmin(1), labels(phlmin(1))])

[ssmax, sslmax] = sort(maxssgain); disp([ssmax(end), labels(sslmax(end))])
[ssmin, sslmin] = sort(minssgain); disp([ssmin(1), labels(sslmin(1))])

disp('***************************************')
% Determine the second maximum and second minimum value and the respective parameters
[trmax, trlmax] = sort(maxtrise); disp([trmax(end-1), labels(trlmax(end-1))])
[trmin, trlmin] = sort(mintrise); disp([trmin(1+1), labels(trlmin(1+1))])

[pwmax, pwlmax] = sort(maxpulsewidth); disp([pwmax(end-1), labels(pwlmax(end-1))])
[pwmin, pwlmin] = sort(minpulsewidth); disp([pwmin(1+1), labels(pwlmin(1+1))])

[phmax, phlmax] = sort(maxpulseheight); disp([phmax(end-1), labels(phlmax(end-1))])
[phmin, phlmin] = sort(minpulseheight); disp([phmin(1+1), labels(phlmin(1+1))])

[ssmax, sslmax] = sort(maxssgain); disp([ssmax(end-1), labels(sslmax(end-1))])
[ssmin, sslmin] = sort(minssgain); disp([ssmin(1+1), labels(sslmin(1+1))])

disp('***************************************')
% Determine the third maximum and third minimum value and the respective parameters
[trmax, trlmax] = sort(maxtrise); disp([trmax(end-2), labels(trlmax(end-2))])
[trmin, trlmin] = sort(mintrise); disp([trmin(1+2), labels(trlmin(1+2))])

[pwmax, pwlmax] = sort(maxpulsewidth); disp([pwmax(end-2), labels(pwlmax(end-2))])
[pwmin, pwlmin] = sort(minpulsewidth); disp([pwmin(1+2), labels(pwlmin(1+2))])

[phmax, phlmax] = sort(maxpulseheight); disp([phmax(end-2), labels(phlmax(end-2))])
[phmin, phlmin] = sort(minpulseheight); disp([phmin(1+2), labels(phlmin(1+2))])

[ssmax, sslmax] = sort(maxssgain); disp([ssmax(end-2), labels(sslmax(end-2))])
[ssmin, sslmin] = sort(minssgain); disp([ssmin(1+2), labels(sslmin(1+2))])

npara = 14;
allpara = [trlmax(npara-2:npara) trlmin(1:3)...
    pwlmax(npara-2:npara) pwlmin(1:3)...
    phlmax(npara-2:npara) phlmin(1:3)...
    sslmax(npara-2:npara) sslmin(1:3)];
for ii = 1:npara
    hy1para(ii) = sum(allpara(:) == ii);
end

angle_para = deg2rad((360/npara):(360/npara):360); % Divide equally no. of parameters by 360 deg

%%% Plotting nominal values
th = linspace(0,2*pi,50);

angles = rad2deg(angle_para);

C = categorical(allpara,1:npara,labels);
figure(553)
h = histogram(C,'barwidth',0.7);
h.FaceColor = [0 0.4470 0.7410];
h.EdgeColor = 'b';
h.LineWidth = 1;
pax = gca;
pax.FontSize = 14;
pax.GridAlpha = 0.3;
ylabel('HY-1')
ylim([0 8])
title('Histogram of influential kinetic parameters')


%% Plotting polar plot

figure(31)
% subplot(2,2,1)
for i = 1:length(maxtrise)
    compasspolar(angle_para(i),maxtrise(i),[0 0.4470 0.7410])
    hold on
    compasspolar(angle_para(i),mintrise(i),[0.4660 0.6740 0.1880])
    hold on
end
hold on
polarplot(th,nominaltrise(1)*ones(1,length(th)),'r','LineWidth',2)
title('Rise time (min)')
pax = gca;
pax.ThetaTick = angles;
pax.ThetaTickLabel = labels;
pax.FontSize = 20;
pax.GridAlpha = 0.3;

figure(32)
for i = 1:length(maxpulsewidth)
    compasspolar(angle_para(i),maxpulsewidth(i),[0 0.4470 0.7410])
    hold on
    compasspolar(angle_para(i),minpulsewidth(i),[0.4660 0.6740 0.1880])
    hold on
end
hold on
polarplot(th,nominalpulsewidth(1)*ones(1,length(th)),'r','LineWidth',2)
title('Pulse Width (min)')
pax = gca;
pax.ThetaTick = angles;
pax.ThetaTickLabel = labels;
pax.FontSize = 20;
pax.GridAlpha = 0.3;


figure(33)
for i = 1:length(maxpulseheight)
    compasspolar(angle_para(i),maxpulseheight(i),[0 0.4470 0.7410])
    hold on
    compasspolar(angle_para(i),minpulseheight(i),[0.4660 0.6740 0.1880])
    hold on
end
hold on
polarplot(th,nominalpulseheight(1)*ones(1,length(th)),'r','LineWidth',2)
title('log_{10}(Pulse Height [pM])')
pax = gca;
pax.ThetaTick = angles;
pax.ThetaTickLabel = labels;
pax.FontSize = 20;
pax.GridAlpha = 0.3;

figure(34)
for i = 1:length(maxssgain)
    compasspolar(angle_para(i),maxssgain(i),[0 0.4470 0.7410])
    hold on
    compasspolar(angle_para(i),minssgain(i),[0.4660 0.6740 0.1880])
    hold on
end
hold on
polarplot(th,nominalssgain(1)*ones(1,length(th)),'r','LineWidth',2)
% rlim([0 15])
title('log_{10}(Steady State Value [fM])')
pax = gca;
pax.ThetaTick = angles;
pax.ThetaTickLabel = labels;
pax.FontSize = 20;
pax.GridAlpha = 0.3;


