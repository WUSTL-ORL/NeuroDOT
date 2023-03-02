function [dim, resampled_A] = resample_A_to_Atlas(A_subj, dim_subj, info_atlas, affineTform)

%   RESAMPLE_A_TO_ATLAS transforms the A matrix of a subject to the atlas
%   space
%   
%   Each measurement in the A matrix is converted to a volume and affine
%   transformed to atlas space. The A matrix is resampled according to the
%   good voxels in the transformed volume. A new info.tissue.dim is
%   generated that matches the resampled A matrix.
%
%   Copyright (c) 2017 Washington University 
%   Created By: Abigail L. Magee
%   Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.

%% Parameters and Initialization.
if ~exist('affineTform', 'var')  ||  isempty(affineTform)
    affineTform = eye(4);
end

channelNo = size(A_subj,1);
maskA = Good_Vox2vol(ones(size(dim_subj.Good_Vox)), dim_subj);
maskA_atlas = affine3d_img(maskA, dim_subj, info_atlas, affineTform, 'nearest');
goodVox = find(maskA_atlas);
resampled_A = zeros(channelNo, length(goodVox));

%% Loop through channels and resample to match atlas
for kk=1:channelNo
    A_meas = A_subj(kk,:);
    vol_meas = Good_Vox2vol(A_meas,dim_subj);
    tformed_vol = affine3d_img(vol_meas, dim_subj, info_atlas, affineTform, 'nearest');
    resampled_A(kk,:) = tformed_vol(goodVox);
end
dim = info_atlas;           % Match dim to atlas
dim.Good_Vox = goodVox;     % Good voxels
dim.nVt = channelNo;        % Number of measurements
dim.sV = info_atlas.mmx;    % Metric size of a voxel
