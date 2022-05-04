function kern=MakeMeKernel(k_size,type)
%
% This function creates a spatial kernel for use with region of interest
% sampling.
% Inputs:
%   k_size  The spatial size of the kernel in voxels. This is the diameter
%           for cube and gaussian but the radius for sphere.
%   type    Option for kernal form. Current options are 'sphere' (default),
%           'cube' and 'gaussian'.
%
% Outputs: 
%   kern    the kernel
%
%
%% Parameters and initialization
if ~exist('k_size','var'),  k_size=5;       end
if ~exist('type','var'),    type='sphere';  end


%% Make kernel
switch type
    case 'gaussian'
        kern = zeros(k_size, k_size, k_size);
        kern(ceil(k_size/2), ceil(k_size/2), ceil(k_size/2)) = 1;
        kern = smooth3(kern, 'gaussian', k_size, 1.2);
    case 'cube'
        kern = ones(k_size, k_size, k_size);
    case 'sphere'
        [xgv, ygv, zgv] = meshgrid(-k_size:k_size, -k_size:k_size, -k_size:k_size);
        kern = sqrt(xgv.^2 + ygv.^2 + zgv.^2) <= k_size;
end