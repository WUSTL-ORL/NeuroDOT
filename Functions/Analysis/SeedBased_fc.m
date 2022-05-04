function [fcMaps,fcMatrix,seedTT]=SeedBased_fc(data,ROI,dim,params)
%
% This function calculates Seed Based functional connectivity between a
% volumetric data set (data) and a set of regions of interest (ROI). The
% ROI can be in 'index' coordinates of the voxels (default) or they can 
% be true coordinates. params.seed_type can be set for either 'idx' or
% 'coord'. dim is the metadata for the volumetric space. To increase
% computational efficiency, the data is assumed to be in 2 dimensions with
% time in the 2nd dimension. If data is not in 2D, it is reshaped to be in
% 2D. The outputs are the Fisher-z transformed maps (fcMaps), the Fisher-z
% transformed fc matrix (fcMatrix), and the ROI timecourses (seedTT). A
% temporal mask delineating time points to keep may also be passed in
% params.tMask.


%% Parameters and initialization
dims = size(data);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data) > 2);
if NDtf
    data = reshape(data, [], Nt);
end

if ~exist('params','var'), params=struct;end
if ~isfield(params,'seed_type'),params.seed_type='idx';end
if ~isfield(params,'seed_rad'),params.seed_rad=5;end
if ~isfield(params,'tMask'),params.tMask=ones(Nt,1);end

if strcmp(params.seed_type,'coord')
    ROI = change_space_coords(ROI, dim, 'idx');
end

kern=floor(params.seed_rad/2);
Nseeds=size(ROI,1);
seedTT=zeros(Nt,Nseeds);
fcMaps=zeros(dim.nVx,dim.nVy,dim.nVz,Nseeds);


%% Make kernel for a sphere
disp(['Forming ROIs'])
[xgv, ygv, zgv] = meshgrid(-kern:kern,-kern:kern,-kern:kern);
kernel = sqrt(xgv.^2 + ygv.^2 + zgv.^2) <= kern;


%% Normalize data 
disp(['Normalizing data'])
DC=squeeze(nanmean(bsxfun(@times,data,params.tMask'),2));
data=bsxfun(@times,data,params.tMask')-repmat(DC,[1,Nt]); % Subtract off mean
data=normrND(data);               % Normalize across time
if NDtf
    data=reshape(data,dims);
else
    data=Good_Vox2vol(data,dim);
end


%% Calculate fc
disp(['Calculating correlations'])
for k=1:Nseeds
    roi=zeros(dim.nVx,dim.nVy,dim.nVz);
    roi(ROI(k,1),ROI(k,2),ROI(k,3))=1;
    roi=convn(roi,kernel,'same');
    Ns=nansum(roi(:));
    
    % Generate Seed Time Trace
    seedTT(:,k)=[sum(reshape(bsxfun(@times,data,roi),[],Nt),1)./Ns]';
    
    % Generate Seed map
    fcMaps(:,:,:,k)=reshape(FisherR2Z((seedTT(:,k)'*...
        reshape(data,[],Nt)')'),size(roi));
end

fcMatrix=FisherR2Z(corrcoef(seedTT));