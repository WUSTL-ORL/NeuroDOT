%% NEURODOT FULL PROCESSING SCRIPT
% This script combines the Preprocessing and Reconstruction pipelines.
% A set of sample Data files exist in the /Data directory of the toolbox.
% For each data set, there is an example *.pptx with selected
% visualizations. The NeuroDOT_Tutorial_Full_Data_Processing.pptx uses
% NeuroDOT_Data_Sample_CCW1.mat. 
%


%% PREPROCESSING PIPELINE

%% Load Measurement data
dataset='CCW1'; % CCW1, CCW2, CW1, IN1, OUT1, GV1, HW1, HW2, HW3_Noisy,  RW1
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
        A_fn='A_Adult_96x92_on_Example_Mesh_test.mat';   % Sensitivity Matrix
        dt=30;                      % Block length
        tp=16;                      % Example (block averaged) time point
        
end


%% PRE-PREOCESSING PIPELINE
lmdata = logmean(data);                                 % Logmean Light Levels
info = FindGoodMeas(lmdata, info, 0.075);               % Detect Noisy Channels
lmdata = detrend_tts(lmdata);                           % Detrend Data
lmdata = highpass(lmdata, .02, info.system.framerate);  % High Pass Filter (0.02 Hz)
lmdata = lowpass(lmdata, 1, info.system.framerate);     % Low Pass Filter 1 (1.0 Hz)
hem = gethem(lmdata, info);                             % Superficial Signal Regression
[lmdata, ~] = regcorr(lmdata, info, hem);
lmdata = lowpass(lmdata, 0.2, info.system.framerate);   % Low Pass Filter 2 (0.5 Hz)
[lmdata, info] = resample_tts(lmdata, info, 1, 1e-5);   % 1 Hz Resampling (1 Hz)


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

%% Block Averaging the measurement data and view
badata = BlockAverage(lmdata, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);

badata=bsxfun(@minus,badata,mean(badata,2));

figure('Position',[100 100 550 780])
subplot(2,1,1); plot(badata(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(badata(keep,:))));
subplot(2,1,2); imagesc(badata(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')


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

% Spectroscopy
if ~exist('E', 'var'),load('E.mat'),end
cortex_Hb = spectroscopy_img(cortex_mu_a, E);
cortex_HbO = cortex_Hb(:, :, 1);
cortex_HbR = cortex_Hb(:, :, 2);
cortex_HbT = cortex_HbO + cortex_HbR;

%% Select Volumetric visualizations of block averaged data
if ~exist('MNI', 'var')
    [MNI,infoB]=LoadVolumetricData('Segmented_MNI152nl_on_MNI111',[],'4dfp'); % load MRI (same data set as in A matrix dim)
end
MNI_dim = affine3d_img(MNI,infoB,A.info.tissue.dim,eye(4),'nearest'); % transform to DOT volume space 

%% Explore the block-averaged data interactively
badata_HbO = BlockAverage(cortex_HbO, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1));
badata_HbOvol = Good_Vox2vol(badata_HbO,A.info.tissue.dim);

Params.Scale=0.8*max(abs(badata_HbOvol(:)));
Params.Th.P=0.2*Params.Scale;
Params.Th.N=-Params.Th.P;
Params.Cmap='jet';
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info)

%% Visualize block-averaged absorption and HbR, HbT data 
% use either of the next two lines to select mu_a (wavelength of 1 or 2 in the third
% dimension) or HbO, HbR, or HbT (second line)
% badata_HbO = BlockAverage(cortex_mu_a(:,:,2), info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata_HbO = BlockAverage(cortex_HbT, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1));
badata_HbOvol = Good_Vox2vol(badata_HbO,A.info.tissue.dim);

Params.Scale=0.8*max(abs(badata_HbOvol(:)));
Params.Th.P=0.2*Params.Scale;
Params.Th.N=-Params.Th.P;
Params.Cmap='jet';
% Params.slices = [49 8 17]; 
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info)

