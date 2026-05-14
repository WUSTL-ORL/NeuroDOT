%% Demonstration of NeuroDOT GLM Functionality

% ------------------------------------------------------------------------
%
% This script will load preprocessed data and run a basic GLM,
% removing high-motion frames that exceed a GVTD threshold 
% and rendering 3 contrast parametric maps on a cortical surface.
%
% The data here is an "Alternating Checkerboard (AC) task -- presentation
% of a black/white checkerboard to either the left or right visual field.
% We expect these to activate the contralaterial visual cortex.
%
% ------------------------------------------------------------------------

% load HRF (convolved with task to create model regressor) and
% cortical mesh (for rendering) -- these are included with NeuroDOT

load('hrf_DOT3.mat');
[MNI,infoMNI]=LoadVolumetricData('Segmented_MNI152nl_on_333',[],'4dfp');

load('MNI164k_big.mat') 
[T1,infoT1] = LoadVolumetricData('mni152nl_T1_on_333', '', '4dfp');

%% Load data

load('AC001_Reconstructed_Data_Sample.mat'); 

% GVTD params

params_AC.GVTD_censor = 1; % 1 = apply GVTD frame deletion
params_AC.GVTD_th= 0.001;  % threshold if doing GVTD frame deletion
params_AC.DoFilter = 1;    % 1 = bandpass filter the data

% run GLM -- see mfile header for input and output description
%
% this will generate a figure displaying i) the design matrix,
% ii) GVTD superimposed on the task time course, and iii) a histogram 
% of GVTD values

[b,e,DM,EDM] = GLM_181206(squeeze(cortex_Hb(:,:,1)),hrf,info,params_AC);

figure('Color','k','Position',[100,100,800,1000]);

% loop over stim type, generate a beta map for each

stimtype = {'left', 'right', 'contrast'};

for j = 1:length(stimtype)
    
    stim = stimtype{j};
    
    switch stim
        case 'left'
            % map of left stimulus
            bL=squeeze(b(:,3)-b(:,2));
            bvol=Good_Vox2vol(bL,info.tissue.dim);
            
        case 'right'
            % map of right stimulus
            bR=squeeze(b(:,4)-b(:,2));
            bvol=Good_Vox2vol(bR,info.tissue.dim);
            
        case 'contrast'
            % map of contrast (left - right)
            bC=squeeze(b(:,3)-b(:,4));
            bvol=Good_Vox2vol(bC,info.tissue.dim);
    end

    % affine transform beta map to MNI space

    map_in_MNI_space = affine3d_img(bvol,info.tissue.dim,infoMNI,info.tissue.affine); 

    % Visualization

    % color scale parameters -- values here give good results

    params_AC.Scale=0.8*max(bvol(bvol(:)~=0));   
    s=params_AC.Scale;
    params_AC.Th.P=0;                       
    params_AC.Th.N=-params_AC.Th.P;
    tp=params_AC.Th.P;
    tn=params_AC.Th.N;

    params_AC.view='post';  % posterior view for visual task
    params_AC.ctx='std';    % standard inflation
    params_AC.OL=0;

    % Surface Render

    subplot(3,1,j);
    params_AC.fig_handle=gca;
    PlotInterpSurfMesh(map_in_MNI_space,MNIl,MNIr,infoMNI,params_AC)
    title(['Beta Map -- AC001/task: ' stimtype{j} ],'Color','w', 'FontSize', 18);

end

% echo rendering params to command window to document

params_AC

