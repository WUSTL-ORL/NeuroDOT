function [mapL,mapR,hPatchesL,hPatchesR,params2] = PlotInterpSurfMesh(volume, meshL, meshR, dim, params)

% PLOTINTERPSURFMESH Interpolates volumetric data onto hemispheric meshes for
% display.
%
%   PLOTINTERPSURFMESH(volume, meshL, meshR, dim) takes functional imaging
%   data "volume" and interpolates it onto the left and right hemispheric
%   surface meshes given in "meshL" and "meshR", respectively, using the
%   spatial information in "dim". The result is a bilateral sagittal view
%   of the activations overlain onto the surface of the brain represented
%   by the meshes.
%
%   PLOTINTERPSURFMESH(volume, meshL, meshR, dim, params) allows the user to
%   specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       Scale       (90% max)       Maximum value to which image is scaled.
%       Th                          Thresholds.
%           .P      (25% max)       Value of min threshold to display
%                                   positive data values.
%           .N      (-Th.P)         Value of max threshold to display
%                                   negative data values.
%
% Dependencies: VOL2SURF_MESH, PLOTLRMESHES, ADJUST_BRAIN_POS, ROTATE_CAP,
% ROTATION_MATRIX.
%
% See Also: PLOTSLICES.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

%% Parameters and Initialization.
if (ndims(volume) >= 4)  &&  (size(volume, 4) > 1)
    % If volume is >=4D, or if so, and the time dimension > 1...
    error(['*** "volume" input has more than one time point.',...
        ' PlotInterpSurfMesh can only image one time point. ***'])
end

if ~exist('params', 'var')  ||  isempty(params)
    params = [];
end


%%  Check Volume Dimensions
[volNx,volNy,volNz] = size(squeeze(volume(:,:,:,1)));

if (volNx < dim.nVx || volNy < dim.nVy || volNz < dim.nVz)
   warning(['*** volume and surface mesh have incompatible dimensions',...
   sprintf('\n\t (volume: %d %d %d ; mesh: %d %d %d)\n',volNx,volNy,volNz,dim.nVx,dim.nVy,dim.nVz),...
   'volume may render incorrectly (or not at all) ***']);
end

if ~exist('params', 'var')  ||  isempty(params)
    params = [];
end


%% Interpolate surface mesh into maps.
mapL = vol2surf_mesh(meshL, volume, dim, params);
mapR = vol2surf_mesh(meshR, volume, dim, params);


%% Image the left and right maps.
[hPatchesL,hPatchesR,params2] = PlotLRMeshes(mapL, mapR, params);



%
