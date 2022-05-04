function tt=ExtractVolumetricROI_TimeTrace(data,Coord,k_size,type)
%
% This function extracts a time trace from 4D data assuming time is the
% last dimension, averaging data within a ROI with a given structure.
% Inputs:
%   data    4D array of data
%   Coord   Nx3 vector containing coordinates for kernel center(s)
%   k_size  The spatial size of the kernel in voxels. This is the diameter
%           for cube and gaussian but the radius for sphere.
%   type    Option for kernal form. Current options are 'sphere' (default),
%           'cube' and 'gaussian'.
%
% Outputs: 
%   tt      The time traces extracted from ROI centered at Coord
%
%
%% Parameters and initialization
if ~exist('k_size','var'),  k_size=5;       end
if ~exist('type','var'),    type='sphere';  end
[Nx,Ny,Nz,Nt]=size(data);
Ntt=size(Coord,1);
tt=zeros(Ntt,Nt);


%% Make kernel
kernVol=MakeVolumetricKernel([Nx,Ny,Nz],Coord,k_size,type);


%% Extract data
for j=1:Ntt
    tempKern=squeeze(kernVol(:,:,:,j));
    dataTT = bsxfun(@times,data,tempKern);
    tt(j,:) = squeeze(sum(sum(sum(dataTT))))./nnz(tempKern(:));
end