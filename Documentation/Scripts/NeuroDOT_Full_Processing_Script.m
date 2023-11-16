%% NEURODOT FULL PROCESSING SCRIPT
% This script combines the Preprocessing and Reconstruction pipelines.
% A set of sample Data files exist in the /Data directory of the toolbox.
% For each data set, there is an example *.pptx with selected
% visualizations. The NeuroDOT_Tutorial_Full_Data_Processing.pptx uses
% NeuroDOT_Data_Sample_CCW1.mat. 

%% Add paths
addpath(genpath('/data/culver/data1/matlab_codes/NIRFASTer-master'),'-BEGIN')
addpath(genpath('/data/culver/data1/matlab_codes/NeuroDOT_NITRC-Release'),'-BEGIN')
addpath(genpath('/data/culver/data1/matlab_codes/NeuroDOT_Internal_Additional_Files_and_Functions'),'-BEGIN')
addpath(genpath('/data/culver/data1/matlab_codes/gifti-1.4'),'-BEGIN')
addpath(genpath('/data/culver/data1/Paul/SCOT_Paul/util'))

%% Load Measurement data
dataset='CCW2'; % CCW1, CCW2, CW1, IN1, OUT1, GV1, HW1, HW2, HW3_Noisy,  RW1
load(['NeuroDOT_Data_Sample_',dataset,'.mat']); % data, info, flags

% Set parameters for A and block length for quick processing examples
switch dataset
    case {'CCW1','CCW2','CW1','OUT'}
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


%% General data QC with synchpts if present
Plot_RawData_Time_Traces_Overview(data,info);           % Time traces
info = Plot_RawData_Cap_DQC(data, info);                % Cap-relevant views
info = Plot_RawData_Metrics_II_DQC(data,info)                  % Raw data quality figs


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


%% PRE-PREOCESSING PIPELINE
% Note: the first 3 lines are repeated from above but with changed variable names
lmdata = logmean(data);                                                   % Logmean Light Levels
info = FindGoodMeas(lmdata, info, 0.075);                                 % Detect Noisy Channels
lmdata = detrend_tts(lmdata);                                             % Detrend Data
lmdata = highpass(lmdata, .02, info.system.framerate);                    % High Pass Filter (0.02 Hz)
lmdata = lowpass(lmdata, omega_lp1, info.system.framerate);                       % Low Pass Filter 1 (1.0 Hz)
hem = gethem(lmdata, info);                                               % Superficial Signal Regression
[lmdata, ~] = regcorr(lmdata, info, hem);
lmdata = lowpass(lmdata, 0.5, info.system.framerate);                     % Low Pass Filter 2 (0.5 Hz)
[lmdata, info] = resample_tts(lmdata, info, 1, 1e-5);                     % 1 Hz Resampling (1 Hz)
[info.GVTD, info.DQ_metrics.med_GVTD] = CalcGVTD(lmdata(info.MEAS.GI & info.pairs.r2d<20,:));         % Calculate GVTD


%% View pre-processed data
keep = info.pairs.WL==2 & info.pairs.r2d < 40 & info.MEAS.GI; % measurements to include

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lmdata(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(lmdata(keep,:))));
subplot(3,1,2); imagesc(lmdata(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')
[ftmag,ftdomain] = fft_tts(squeeze(mean(lmdata(keep,:),1)),info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag);
xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])

nlrGrayPlots_180818(lmdata,info); % Gray Plot with synch points

Plot_TimeTrace_With_PowerSpectrum(lmdata,info); % As above, but now automated with all wavelengths
GrayPlots_Rsd_by_Wavelength(lmdata,info); % As above but now with all wavelengths


%% Block Averaging the measurement data and view
badata = BlockAverage(lmdata, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata=bsxfun(@minus,badata,mean(badata(:,1:4),2));

figure('Position',[100 100 550 780])
subplot(2,1,1); plot(badata(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(badata(keep,:))));
subplot(2,1,2); imagesc(badata(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');
ylabel('Measurement #')


%% RECONSTRUCTION PIPELINE
if ~exist('A', 'var')       % In case running by hand or re-running script
    A=load([A_fn],'info','A');
    if length(size(A.A))>2  % A data structure [wl X meas X vox]-->[meas X vox]
        [Nwl,Nmeas,Nvox]=size(A.A);
        A.A=reshape(permute(A.A,[2,1,3]),Nwl*Nmeas,Nvox);
    end        
end
Nvox=size(A.A,2);
Nt=size(lmdata,2);
cortex_mu_a=zeros(Nvox,Nt,2);
for j = 1:2
    keep = (info.pairs.WL == j) & (info.pairs.r2d <= 40) & info.MEAS.GI;
    disp('> Inverting A')                
    iA = Tikhonov_invert_Amat(A.A(keep, :), 0.01, 0.1); % Invert A-Matrix
    disp('> Smoothing iA')
    iA = smooth_Amat(iA, A.info.tissue.dim, 3);         % Smooth Inverted A-Matrix      
    cortex_mu_a(:, :, j) = reconstruct_img(lmdata(keep, :), iA);% Reconstruct Image Volume
end


%% Spectroscopy
if ~exist('E', 'var'),load('E.mat'),end
cortex_Hb = spectroscopy_img(cortex_mu_a, E);
cortex_HbO = cortex_Hb(:, :, 1);
cortex_HbR = cortex_Hb(:, :, 2);
cortex_HbT = cortex_HbO + cortex_HbR;


%% Select Volumetric visualizations of block averaged data
if ~exist('MNI', 'var')
[MNI,infoB]=LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti',[],'nii'); % Load MRI (same data set as in A matrix dim)
end
MNI_dim = affine3d_img(MNI,infoB,A.info.tissue.dim,eye(4),'nearest'); % Transform to DOT volume space

% Block Average Data
badata_HbO = BlockAverage(cortex_HbO, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1));
badata_HbOvol = Good_Vox2vol(badata_HbO,A.info.tissue.dim);
tp_Eg=squeeze(badata_HbOvol(:,:,:,tp));

% Explore PlotSlices - The basics (Slide 22 in ppt)
PlotSlices(MNI_dim)                             % Anatomy only
PlotSlices(MNI_dim,A.info.tissue.dim)           % Anatomy + volumetric data
PlotSlices(MNI_dim,A.info.tissue.dim,params,tp_Eg); % Anatomy + volumetric data + functional data

% Visualize the data (Slide 23 in ppt)
PlotSlices(tp_Eg,A.info.tissue.dim);            % Data by itself
PlotSlices(MNI_dim,A.info.tissue.dim,[],tp_Eg); % Data with anatomical underlay
% Set parameters to visualize more specific aspects of data
Params.Scale=0.8*max(abs(tp_Eg(:)));     % Scale wrt/max of data
Params.Th.P=0.4*Params.Scale;            % Threshold to see strong activations
Params.Th.N=-0.010;                % Thresholds go both ways
Params.Cmap='jet';
PlotSlices(MNI_dim,A.info.tissue.dim,Params,tp_Eg);

% Explore the block-averaged data a bit more interactively (slide 24 in ppt)
Params.Scale=0.8*max(abs(badata_HbOvol(:)));
Params.Th.P=0;
Params.Th.N=-Params.Th.P;
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info)

% Explore the not-block-averaged data a bit more interactively (slide 24 in ppt)
HbOvol = Good_Vox2vol(cortex_HbO,A.info.tissue.dim);
Params.Scale=1e-3;
Params.Th.P=1e-4;
Params.Th.N=-Params.Th.P;
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,HbOvol,info)


%% Select Surface visualizations
if ~exist('MNIl', 'var'),load(['MNI164k_big.mat']);end
HbO_atlas = affine3d_img(badata_HbOvol,A.info.tissue.dim,infoB,eye(4));
tp_Eg_atlas=squeeze(HbO_atlas(:,:,:,tp));
pS=Params;
pS.view='post';

pS.ctx='std'; % Standard pial cortical view
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

pS.ctx='inf'; % Inflated pial cortical view
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

pS.ctx='vinf';% Very Inflated pial cortical view
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

