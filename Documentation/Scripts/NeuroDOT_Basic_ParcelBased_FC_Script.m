
%% NeuroDOT Functional Connectivity Example Script

% ------------------------------------------------------------------------
%
% This script will demonstrate generation of functional connectivity
% maps using the Gordon cortical parcellation. It uses prepared data:
% Subj3_Nph2014_REST-image.mat, which must be downloaded and placed
% in an accessible directory (the script will change the Matlab working
% directory to this location). 
%
% ------------------------------------------------------------------------

%% Path Setup
%  set f0 to the path where Subj3_Nph2014_REST-image.mat lives

% % f0 = '/Path/to/Data/';
% % cd(f0);

%% Load Data
%  this will create variables "cortex_Hb" and "info" in workspace

load('Subj3_Nph2014_REST-image.mat');  

%% Load Support Files

% these files are included with NeuroDOT

[ MNI_parcels,parcelsInfo ] = LoadVolumetricData('Parcels_MNI_333','','4dfp');

load('MNI_coord_meshes_32k', 'MNIl', 'MNIr');
load('IM_Gordon_2014_333_Parcels.mat','IM','Parcels','Parcel_Nets');
load('MNI164k_big.mat');

% these files are not included with NeuroDOT and require separate download

load('Adult_96x92_Overlapped_Gordon_Parcels.mat');
load('Adult_96x92_FOV.mat');
A = load('A_Adult_96x92.mat','info');

%% Prepare HbO data for FC calculation

% set up GVTD frame censoring

GVTDth = 1e-3; 
tMask = info.GVTD<GVTDth;
Nt2 = sum(tMask);

% prepare data

HbO = cortex_Hb(:,:,1);                 
HbOtm = HbO(:,tMask);                         % exclude high motion frames
HbOm = squeeze(mean(HbOtm,2));                % find DC (global mean)
HbOv = HbOtm-repmat(HbOm,[1,Nt2]);            % Subtract off mean
HbOv = normr(HbOv);                           % Normalize across time
HbOv = Good_Vox2vol(HbOv,A.info.tissue.dim);  % convert nvoxels x time to (nVx x nVy x nVz) x time

% put data in MNI space to match template
HbOvol=affine3d_img(HbOv,A.info.tissue.dim,parcelsInfo,A.info.tissue.affine);

%% Calculate Parcel-based Functional Connectivity
% this function may take several minutes to run

[ fcMaps,fcMatrix,ParcelTT ] = ParcelBased_fc(HbOvol,Overlap_parcel,parcelsInfo);
fcMaps = fcMaps.*FOV;

%% Plot example FC maps for parcels of interest
%  Visual, Auditory, SMhand, and Fronto-Parietal Parcels

% because only parcels that intersect the array FOV are included in the
% analysis, the parcel numbering here differs from the original assignments
% in the Gordon atlas. See ParcelOrder.xlxs for a complete listing.

visualParcels       = [42,75];
auditoryParcels     = [24,96];
somatomotorParcels  = [11,86];
fpParcels           = [51,77];

parcelList = [visualParcels; auditoryParcels; somatomotorParcels; fpParcels];

figure('Color','k','Position',[100,100,1600,500]);

for i = 1:8
    subplot(2,4,i)
    parcelNum = parcelList(i);
    fcMap = fcMaps(:,:,:,parcelNum);
    params.Scale = 0.5*max(fcMap(:));
    params.Th.P=0;params.Th.N=0;
    title(['Parcel ', num2str(parcelNum)],'Color','w');
    params.fig_handle = gca;
    PlotInterpSurfMesh(squeeze(fcMap),MNIl,MNIr,parcelsInfo,params);
end

ha = axes('Position',[0 0 1 1],'Visible','off');
text(0.2, 0.98, 'Visual','Units','normalized','HorizontalAlignment','left', ...
     'FontWeight','bold','Color','w','Parent',ha);
text(0.4, 0.98, 'Auditory','Units','normalized','HorizontalAlignment','left', ...
 'FontWeight','bold','Color','w','Parent',ha);
text(0.6, 0.98, 'SMhand','Units','normalized','HorizontalAlignment','left', ...
 'FontWeight','bold','Color','w','Parent',ha);
text(0.8, 0.98, 'Fronto-Parietal','Units','normalized','HorizontalAlignment','left', ...
 'FontWeight','bold','Color','w','Parent',ha);

