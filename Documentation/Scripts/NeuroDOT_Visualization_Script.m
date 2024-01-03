edit%% NEURODOT VISUALIZATION TUTORIAL SCRIPT
% This script uses the NeuroDOT_Full_Processing_Script as a base to help
% teach new users about the different visualization options available in
% the NeuroDOT toolbox
% A set of sample Data files exist in the /Data directory of the toolbox.
% For each data set, there is an example *.pptx with selected
% visualizations. The NeuroDOT_Tutorial_Visualization.pptx uses
% NeuroDOT_Data_Sample_CCW1.mat. 



%% Load Measurement data
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


%% General data QC with synchpts if present
Plot_RawData_Time_Traces_Overview(data,info);           % Time traces
Plot_RawData_Cap_DQC(data,info);                        % Cap-relevant views
Plot_RawData_Metrics_II_DQC(data,info)                  % Raw data quality figs


%% Pre-process data
% Pre-processing pipeline, but only the part before filtering data
lmdata_b4_filt = logmean(data);                                           % Logmean Light Levels
info_b4_filt  = FindGoodMeas(lmdata_b4_filt , info, 0.075);               % Detect Noisy Channels
lmdata_b4_filt  = detrend_tts(lmdata_b4_filt );                           % Detrend Data
info_b4_filt.GVTD = CalcGVTD(lmdata_b4_filt(info_b4_filt.MEAS.GI & info_b4_filt.pairs.r2d<20,:));         % Calculate GVTD

% Full Pre-processing Pipeline
lmdata = logmean(data);                                                   % Logmean Light Levels
info = FindGoodMeas(lmdata, info, 0.075);                                 % Detect Noisy Channels
lmdata = detrend_tts(lmdata);                                             % Detrend Data
lmdata = highpass(lmdata, .02, info.system.framerate);                    % High Pass Filter (0.02 Hz)
lmdata = lowpass(lmdata, 1, info.system.framerate);                       % Low Pass Filter 1 (1.0 Hz)
hem = gethem(lmdata, info);                                               % Superficial Signal Regression
[lmdata, ~] = regcorr(lmdata, info, hem);
lmdata = lowpass(lmdata, 0.5, info.system.framerate);                     % Low Pass Filter 2 (0.5 Hz)
[lmdata, info] = resample_tts(lmdata, info, 1, 1e-5);                     % 1 Hz Resampling (1 Hz)
info.GVTD = CalcGVTD(lmdata(info.MEAS.GI & info.pairs.r2d<20,:));         % Calculate GVTD


%% Visualize Effects of Pre-Proecssing
% Visualize Pre-filtering
% Measurements to include
keep_b4_filt = info_b4_filt.pairs.WL==2 & info_b4_filt.pairs.r2d < 40 & info_b4_filt.MEAS.GI; 

% Visualize
figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lmdata_b4_filt(keep_b4_filt,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m_b4_filt=max(max(abs(lmdata_b4_filt(keep_b4_filt,:))));
subplot(3,1,2); imagesc(lmdata_b4_filt(keep_b4_filt,:),[-1,1].*m_b4_filt); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')
[ftmag_b4_filt,ftdomain_b4_filt] = fft_tts(squeeze(mean(lmdata_b4_filt(keep_b4_filt,:),1)),info_b4_filt.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain_b4_filt,ftmag_b4_filt);
xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])

nlrGrayPlots_180818(lmdata_b4_filt,info_b4_filt); % Gray Plot with synch points


% Visualize filtered data
% Measurements to include
keep = info.pairs.WL==2 & info.pairs.r2d < 40 & info.MEAS.GI; % measurements to include

% Visualize
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


% Compare 750nm and 850nm wavelengths
Plot_TimeTrace_With_PowerSpectrum(lmdata,info); % Data, Grayplots organized by measurement #, Powerspectrum
GrayPlots_Rsd_by_Wavelength(lmdata,info); %Grayplots organized by SD separation


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

% Explore PlotSlices - The basics
PlotSlices(MNI_dim)                             % Anatomy only
PlotSlices(MNI_dim,A.info.tissue.dim)           % Anatomy + volumetric data
PlotSlices(MNI_dim,A.info.tissue.dim,[],tp_Eg); % Anatomy + volumetric data + functional data

% Match your view of PlotSlices to the images in the tutorial PPT (slide 27)
PlotSlices(MNI_dim,A.info.tissue.dim,[],tp_Eg);

% Explore PlotSlices - in depth (Slides 28-31)
% Defaults only, no params structure as input
PlotSlices(MNI_dim,A.info.tissue.dim,[],tp_Eg);
% Params 1
Params1.PD = 1; % Turn Positive Definite On
PlotSlices(MNI_dim,A.info.tissue.dim,Params1,tp_Eg);
% Params 2
Params2.Th.P = 0.001;
Params2.Th.N = - Params2.Th.P;
PlotSlices(MNI_dim,A.info.tissue.dim,Params2,tp_Eg);
% Params 3
Params3.Th.P = 0.001;
Params3.Th.N = -0.004;
PlotSlices(MNI_dim,A.info.tissue.dim,Params3,tp_Eg);
% Params 4
Params4.Scale = 2*max(abs(tp_Eg(:)));
PlotSlices(MNI_dim,A.info.tissue.dim,Params4,tp_Eg);
% Params 5
Params5.Cmap = 'cool';
PlotSlices(MNI_dim,A.info.tissue.dim,Params5,tp_Eg);
% Params 6
Params6.Cmap.P = 'winter';
Params6.Cmap.N = 'summer';
PlotSlices(MNI_dim,A.info.tissue.dim,Params6,tp_Eg);
% Params 7
Params7.cbmode = 1;
Params7.cbticks = [0,0.25,0.5,0.75,1];
PlotSlices(MNI_dim,A.info.tissue.dim,Params7,tp_Eg);
% Params 8
Params8.cbmode = 1;
Params8.cblabels = {'Hello', 'Goodbye'};
PlotSlices(MNI_dim,A.info.tissue.dim,Params8,tp_Eg);
% Params 9 
Params9.cbmode = 1;
Params9.cbticks = [0,0.25,0.5,0.75,1];
Params9.cblabels = {'1', '2', '3', '4', '5'};
Params9.Cmap.P = 'winter';
Params9.Cmap.N = 'summer';
Params9.Scale = 0.8*max(abs(tp_Eg(:)));
Params9.Th.P = 0.001;
Params9.Th.N = -0.004;
PlotSlices(MNI_dim,A.info.tissue.dim,Params9,tp_Eg);


% Visualize the data (slide 32)
PlotSlices(tp_Eg,A.info.tissue.dim);            % Data by itself
PlotSlices(MNI_dim,A.info.tissue.dim,[],tp_Eg); % Data with anatomical underlay
% Set parameters to visualize more specific aspects of data
Params.Scale=0.8*max(abs(tp_Eg(:)));     % Scale wrt/max of data
Params.Th.P=0.5*Params.Scale;            % Threshold to see strong activations
Params.Th.N=-Params.Th.P;                % Thresholds go both ways
Params.Cmap='jet';
PlotSlices(MNI_dim,A.info.tissue.dim,Params,tp_Eg);

% Explore the block-averaged data a bit more interactively
Params.Scale=0.8*max(abs(badata_HbOvol(:)));
Params.Th.P=0;
Params.Th.N=-Params.Th.P;
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info)

% Explore the not-block-averaged data a bit more interactively
HbOvol = Good_Vox2vol(cortex_HbO,A.info.tissue.dim);
Params.Scale=4e-3;
Params.Th.P=1e-3;
Params.Th.N=-Params.Th.P;
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,HbOvol,info)


%% Select Surface visualizations
if ~exist('MNIl', 'var'),load(['MNI164k_big.mat']);end
HbO_atlas = affine3d_img(badata_HbOvol,A.info.tissue.dim,infoB,eye(4));
tp_Eg_atlas=squeeze(HbO_atlas(:,:,:,tp));
pS=Params;
pS.view='post'; % Posterior view
pS.ctx='std'; % Standard pial cortical view
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS); %(Slide 34)


% Play around with visualization parameters
% View (Slide 35)
% Dorsal view
pS.view = 'dorsal';
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Lateral view
pS.view = 'lat';
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

% Inflation (Slide 35)
pS.view='post'; %reset to posterior view before visualizing
% Standard pial cortical view
pS.ctx='std'; 
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Inflated pial cortical view
pS.ctx='inf'; 
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Very Inflated pial cortical view
pS.ctx='vinf';
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

% Thresholds (Slide 36)
pS.ctx='std'; %reset to standard inflation before visualizing
% Set to 0
pS.Th.P = 0;
pS.Th.N = -pS.Th.P;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Set a bit higher
pS.Th.P = 1e-4;
pS.Th.N = -pS.Th.P;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Set even higher
pS.Th.P = 1e-3;
pS.Th.N = -pS.Th.P;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Set very high
pS.Th.P = 1e-2;
pS.Th.N = -pS.Th.P;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

% Scale (Slide 37)
pS.Th.P = 0;
pS.Th.N = -pS.Th.P; %reset thresholds first
% Set to 4e-3
pS.Scale = 4e-3;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Set lower
pS.Scale = 1e-3;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Set higher
pS.Scale = 9e-3;
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);

% Colormaps (Slide 38)
% reset scale and thresholds first
pS.Scale = 4e-3; 
pS.Th.P = 1e-3;
pS.Th.N = -pS.Th.P;
% Same colormap for P/N
pS.Cmap = 'cool';
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);
% Differnt colormaps for P/N
pS = rmfield(pS, 'Cmap'); %remove Cmap field so we can replace it with a struct containing both positive and negative colormaps
pS.Cmap.P = 'winter';
pS.Cmap.N = 'summer';
PlotInterpSurfMesh(tp_Eg_atlas, MNIl,MNIr, infoB, pS);


%% True Color Parameter
params.Cmap = 'jet';
% Get data at another timepoint
tp_Eg_atlas_2=squeeze(HbO_atlas(:,:,:,tp+9)); %adding 9 to tp gets us data when the checkerboard is in a different quadrant

% Make mask of timepoint = 16 data
mask_temp_1 = tp_Eg_atlas;
mask1 = zeros(size(mask_temp_1));
mask1(mask_temp_1 > pS.Th.P) = 1; %set tp1 positive activations to 1
PlotInterpSurfMesh(mask1, MNIl,MNIr, infoB, pS);
mask_temp_2 = tp_Eg_atlas_2;
mask2 = zeros(size(mask_temp_2));
mask2(mask_temp_2 > pS.Th.P) = 1; %set tp2 positive activations to 1
PlotInterpSurfMesh(mask2, MNIl,MNIr, infoB, pS);

% Combine masks to plot
mask3 = zeros(size(mask_temp_1));
mask3(mask1 == 1 & mask2 == 1) = 3;
mask3(mask1 == 1 & mask3 ~= 3) = 1;
mask3(mask2 == 1 & mask3 ~= 3) = 2;

% Set parameters
pS_TC = pS;
pS_TC.Cmap = [1,0,0; 0,0,1; 1,0,1];
pS_TC.TC = 1; %True Color On
pS_TC.OL = 1; %Overlap On
PlotInterpSurfMesh(mask3, MNIl,MNIr, infoB, pS_TC); 





