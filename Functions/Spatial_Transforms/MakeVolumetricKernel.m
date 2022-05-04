function kernVol=MakeVolumetricKernel(VolSize,Coord,k_size,type)
%
% This function makes a kernel of a given types and size and places it in
% the volume at a specific coordinate.
% 
% Inputs:
%   VolSize 1x3 vector containing voxel numbers per dimension
%   Coord   1x3 vector containing coordinates for kernel center(s)
%   k_size  The spatial size of the kernel in voxels. This is the diameter
%           for cube and gaussian but the radius for sphere.
%   type    Option for kernal form. Current options are 'sphere' (default),
%           'cube' and 'gaussian'.
%
% Outputs: 
%   kernVol The volumetric kernel
%
%
%% Parameters and initialization
if ~exist('k_size','var'),  k_size=5;       end
if ~exist('type','var'),    type='sphere';  end
Nkern=size(Coord,1);
kernVol=zeros(VolSize(1),VolSize(2),VolSize(3),Nkern);

%% Make kernel
kern=MakeMeKernel(k_size,type);


%% Place in the volume at the coordinate location
for j=1:Nkern
    kernel = zeros(VolSize);
    kernel(Coord(1), Coord(2), Coord(3)) = 1;
    kernVol(:,:,:,j) = convn(kernel, kern, 'same');
end