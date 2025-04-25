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
dataset='CCW1'; % CCW1, CCW2, CW1, IN1, OUT1, GV1, HW1, HW2, HW3_Noisy,  RW1
load(['NeuroDOT_Data_Sample_',dataset,'.mat']); % data, info, flags

% Set parameters for A and block length for quick processing examples
switch dataset
    case {'CCW1','CCW2','CW1','OUT1'}
        A_fn='A_AdultV24x28.mat';   % Sensitivity Matrix
        dt=36;                      % Block length
        tp=16;                      % Example (block averaged) time point
        
    case {'IN1'}
        A_fn='A_AdultV24x28.mat';   % Sensitivity Matrix
        dt=36;                      % Block length
        tp=32;                      % Example (block averaged) time point
        
    case {'HW1','HW2','RW1','GV1','HW3_Noisy'}
        A_fn='A_Adult_96x92.mat';   % Sensitivity Matrix
        dt=30;                      % Block length
        tp=16;                      % Example (block averaged) time point
        
end



%% General data Quality Assessment with synchpts if present
Plot_RawData_Time_Traces_Overview(data,info);   % Time traces
info = Plot_RawData_Metrics_II_DQC(data,info);         % Spectrum, falloff, and good signal metric
info = Plot_RawData_Cap_DQC(data,info);                % Cap-relevant views


%% View data before filtering
omega_lp1 = 1;

if info.system.framerate/2 < omega_lp1 
    % Adjust Lowpass filter cutoff
    % frequency for systems with lower framerates
    omega_lp1 = (info.system.framerate/2)*0.90;
end

lmdata_b4_filt = logmean(data);                                           % Logmean Light Levels
info_b4_filt  = FindGoodMeas(lmdata_b4_filt , info, 0.075);               % Detect Noisy Channels
lmdata_b4_filt  = detrend_tts(lmdata_b4_filt );                           % Detrend Data
info_b4_filt.GVTD = CalcGVTD(lmdata_b4_filt(info_b4_filt.MEAS.GI & info_b4_filt.pairs.r2d<20,:));         % Calculate GVTD

% Measurements to include
keep = info_b4_filt.pairs.WL==2 & info_b4_filt.pairs.r2d < 40 & info_b4_filt.MEAS.GI; 

% Visualize
figure('Position',[100 100 550 780])
subplot(3,1,1); semilogy(data(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('\phi') 
m=max(max(abs(lmdata_b4_filt(keep,:))));
subplot(3,1,2); imagesc(lmdata_b4_filt(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')
[ftmag,ftdomain] = fft_tts(squeeze(mean(lmdata_b4_filt(keep,:),1)),info_b4_filt.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag);
xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])

nlrGrayPlots_180818(lmdata_b4_filt,info_b4_filt); % Gray Plot with synch points


%% Logmean Light Levels
[lmdata, info.MEAS.Phi_0] = logmean(data);

%% Detect Noisy Channels
info = FindGoodMeas(lmdata, info, 0.075);

% Example visualization
keep = info.pairs.WL==2 & info.pairs.r2d < 40 & info.MEAS.GI; % measurements to include

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lmdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
m=max(max(abs(lmdata(keep,:))));
subplot(3,1,2); imagesc(lmdata(keep,:),[-1,1].*m), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lmdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%% Show nn1, nn2, nn3 (plots)
keepd1=info.MEAS.GI & info.pairs.r2d<20 & info.pairs.WL==2;
keepd2=info.MEAS.GI & info.pairs.r2d>=20 & info.pairs.r2d<30 & info.pairs.WL==2;
keepd3=info.MEAS.GI & info.pairs.r2d>=30 & info.pairs.r2d<40 & info.pairs.WL==2;

figure('Position',[100 100 1000 700])
subplot(3,1,1), plot(lmdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(lmdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(lmdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([500 600]), end

%% Detrend and High-pass Filter the Data
ddata = detrend_tts(lmdata);

% High Pass Filter
hpdata = highpass(ddata, 0.02, info.system.framerate);
% hpdata = highpass(ddata, 0.05, info.system.framerate); % problematic cutoff frequency example

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(hpdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
m=max(max(abs(hpdata(keep,:))));
subplot(3,1,2); imagesc(hpdata(keep,:),[-1,1].*m), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(hpdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[100 100 1000 700])
subplot(3,1,1), plot(hpdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(hpdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(hpdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([500 600]), end


%% Low Pass Filter 1

lp1data = lowpass(hpdata, 1, info.system.framerate);

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lp1data(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
m=max(max(abs(lp1data(keep,:))));
subplot(3,1,2); imagesc(lp1data(keep,:),[-1,1].*m), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lp1data(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[100 100 1000 700])
subplot(3,1,1), plot(lp1data(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(lp1data(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(lp1data(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([500 600]), end
pause
for i = 1:3, subplot(3,1,i), axis auto, end


%% Superficial Signal Regression
hem = gethem(lp1data, info);
[SSRdata, ~] = regcorr(lp1data, info, hem);
% SSRdata = lp1data; % example to ignore SSR

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(SSRdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
m=max(max(abs(SSRdata(keep,:))));
subplot(3,1,2); imagesc(SSRdata(keep,:),[-1,1].*m), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(SSRdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])

%%
figure('Position',[100 100 800 430])
subplot(2,1,1)
plot(hem(2,:))
title('Estimated common superficial signal')
xlabel('Time (samples)')
subplot(2,1,2)
[ftmag,ftdomain] = fft_tts(hem(2,:),info.system.framerate); 
semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])


%%
figure('Position',[100 100 1000 700])
subplot(3,1,1), plot(SSRdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(SSRdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(SSRdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([500 600]), end

%% Low Pass Filter 2
lp2data = lowpass(SSRdata, 0.5, info.system.framerate);
% lp2data = lowpass(SSRdata, 0.05, 10); % example

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lp2data(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(lp2data(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(lp2data(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 10])


%% 1 Hz Resampling
[rdata, info] = resample_tts(lp2data, info, 1, 1e-5);

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(rdata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(3,1,2); imagesc(rdata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image
[ftmag,ftdomain] = fft_tts(mean(rdata(keep,:)',2)',info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag), xlabel('Frequency (Hz)'), ylabel('|X(f)|') % plot vs. log frequency
xlim([1e-3 1])

%%
figure('Position',[100 100 450 700])
subplot(3,1,1), plot(rdata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(rdata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(rdata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')
pause
for i = 1:3, subplot(3,1,i), xlim([100 200]), end


%% Global variance of the temporal derivative (GVTD)
[info.GVTD, info.DQ_metrics.med_GVTD] = CalcGVTD(lp2data(info.MEAS.GI & info.pairs.r2d<20,:));         % Calculate GVTD
info.GVTD_filt_rs=resample_tts(info.GVTD',info,1,1e-5,...
    info.system.init_framerate)';
nlrGrayPlots_180818(rdata,info)

%% Block Averaging
badata = BlockAverage(rdata, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata=bsxfun(@minus,badata,mean(badata(:,1:4),2));
preprocessed = badata;

figure('Position',[100 100 550 780])
subplot(2,1,1); plot(badata(keep,:)'), set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), ylabel('log(\phi/\phi_0)') % plot signals 
subplot(2,1,2); imagesc(badata(keep,:)), colorbar('Location','northoutside'), xlabel('Time (samples)'), ylabel('Measurement #') % show signals as image


%%
figure('Position',[100 100 300 700])
subplot(3,1,1), plot(badata(keepd1,:)'), ylabel('R_{sd} < 20 mm')
subplot(3,1,2), plot(badata(keepd2,:)'), ylabel('R_{sd} \in [20 30] mm')
subplot(3,1,3), plot(badata(keepd3,:)'), ylabel('R_{sd} \in [30 40] mm'), xlabel('Time (samples)')


