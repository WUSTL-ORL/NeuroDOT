%% NEURODOT PREPROCESSING SCRIPT
% This script includes details on the Preprocessing pipeline.
% A file of sample data is already designated below, but you can use the
% "load" command to load your own optical data. In order to load the sample
% file, change the path below in the "addpath" line to the folder under
% which you have ND2 installed.

%% Installation
% installpath = ''; % INSERT YOUR DESIRED INSTALL PATH HERE
% addpath(genpath(installpath))

%% PREPROCESSING PIPELINE


%% S-D Measurements
clear; close all; clc
load('NeuroDOT_Data_Sample_CCW1.mat'); % data, info, flags


%% General data Quality Assessment with synchpts if present
Plot_RawData_Time_Traces_Overview(data,info);   % Time traces
Plot_RawData_Metrics_II_DQC(data,info);         % Spectrum, falloff, and good signal metric
Plot_RawData_Cap_DQC(data,info);                % Cap-relevant views

%% Logmean Light Levels
lmdata = logmean(data);

%% Detect Noisy Channels
info = FindGoodMeas(lmdata, info, 0.075);

% Example visualization
keep = info.pairs.WL==2 & info.pairs.r2d < 40 & info.MEAS.GI; % measurements to include

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(lmdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(lmdata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lmdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%% Show nn1, nn2, nn3 (plots)
info.GVTD = CalcGVTD(lmdata);
nlrGrayPlots_180818(lmdata,info)

keepd1=info.MEAS.GI & info.pairs.r2d<20 & info.pairs.WL==2;
keepd2=info.MEAS.GI & info.pairs.r2d>=20 & info.pairs.r2d<30 & info.pairs.WL==2;
keepd3=info.MEAS.GI & info.pairs.r2d>=30 & info.pairs.r2d<40 & info.pairs.WL==2;

figure('Position',[700 300 1000 700])
subplot(3,1,1), plot(lmdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(lmdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(lmdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([2500 2600]), end

%% Detrend and High-pass Filter the Data
ddata = detrend_tts(lmdata);

% High Pass Filter
hpdata = highpass(ddata, 0.02, info.system.framerate);
% hpdata = highpass(ddata, 0.05, info.system.framerate); % problematic cutoff frequency example

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(hpdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(hpdata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(hpdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[700 300 1000 700])
subplot(3,1,1), plot(hpdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(hpdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(hpdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([2500 2600]), end


%% Low Pass Filter 1
lp1data = lowpass(hpdata, 1, info.system.framerate);

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(lp1data(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(lp1data(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lp1data(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[700 300 1000 700])
subplot(3,1,1), plot(lp1data(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(lp1data(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(lp1data(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([2500 2600]), end
pause
for i = 1:3, subplot(3,1,i), axis auto, end


%% Superficial Signal Regression
hem = gethem(lp1data, info);
[SSRdata, ~] = regcorr(lp1data, info, hem);
% SSRdata = lp1data; % example to ignore SSR

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(SSRdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(SSRdata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(SSRdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[700 600 800 430])
subplot(2,1,1)
plot(hem(2,:))
title('Estimated common superficial signal')
xlabel('Time (samples)')
subplot(2,1,2)
[ftmag,ftdomain] = fft_tts(hem(2,:),info.system.framerate); 
semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])


%%
figure('Position',[700 300 1000 700])
subplot(3,1,1), plot(SSRdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(SSRdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(SSRdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([1600 1600+360*3]), end

%% Low Pass Filter 2
lp2data = lowpass(SSRdata, 0.5, info.system.framerate);
% lp2data = lowpass(SSRdata, 0.05, 10); % example

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(lp2data(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(lp2data(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lp2data(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])


%% 1 Hz Resampling
[rdata, info] = resample_tts(lp2data, info, 1, 1e-5);

figure('Position',[300 600 550 780])
subplot(3,1,1); plot(rdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(rdata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(rdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 1])

%%
figure('Position',[700 300 450 700])
subplot(3,1,1), plot(rdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(rdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(rdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([160 160+36*3]), end


%% Block Averaging
badata = BlockAverage(rdata, info.paradigm.synchpts(info.paradigm.Pulse_2), 36);

preprocessed = badata;

figure('Position',[300 600 550 780])
subplot(2,1,1); plot(badata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(2,1,2); imagesc(badata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
% [ftmag,ftdomain] = fft_tts(mean(badata(keep,:)',2)',info.system.framerate); % Generate average spectrum
% subplot(1,3,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
% xlim([1e-3 10])

%%
figure('Position',[700 300 300 700])
subplot(3,1,1), plot(badata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(badata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(badata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')


