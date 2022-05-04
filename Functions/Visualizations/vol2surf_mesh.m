function Smesh = vol2surf_mesh(Smesh, volume, dim, params)

% VOL2SURF_MESH Interpolates volumetric data onto a surface mesh.
% 
%   Smesh = VOL2SURF_MESH(mesh_in, volume, dim) takes the mesh "Smesh"
%   and interpolates the values of the volumetric data "volume" at the
%   mesh's surface, using the spatial information in "dim". These values
%   are output as "Smesh".
% 
%   Smesh = VOL2SURF_MESH(Smesh, volume, dim, params) allows the user
%   to specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       OL      0   If "overlap" data is presented (OL==1), this sets the
%                   interpolation method to "nearest". Default is "linear".
% 
% See Also: PLOTINTERPSURFMESH, GOOD_VOX2VOL, AFFINE3D_IMG.
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
Ncols = size(volume, 4);
Ncoords = size(Smesh.nodes, 1);
Smesh.data = zeros(Ncoords, Ncols);
extrapval = 0;

nVx = dim.nVx;
nVy = dim.nVy;
nVz = dim.nVz;
dr = dim.mmppix;
center = dim.center;

if ~exist('params', 'var')  ||  isempty(params)
    params = [];
end

if ~isfield(params, 'OL')  ||  isempty(params.OL)
    params.OL = 0;
end
if params.OL
    method = 'nearest';
else
    method = 'linear';
end

buffer = 2;


%% Define coordinate space of volumetric data
X = (dr(1).*[nVx:-1:1]-center(1))';
Y = (dr(2).*[nVy:-1:1]-center(2))';
Z = (dr(3).*[nVz:-1:1]-center(3))';


%% Get coordinates for surface mesh
x = Smesh.nodes(:, 1);
y = Smesh.nodes(:, 2);
z = Smesh.nodes(:, 3);


%% Correct for nodes just otuside of volume (MNI and TT atlas space cuts off occipital pole and part of dorsal tip and lateral extremes).
x((x < (min(X))) & (x > (min(X)-buffer))) = min(X); % 20220317
x((x > (max(X))) & (x < (max(X)+buffer))) = max(X);

y((y < (min(Y))) & (y > (min(Y)-buffer))) = min(Y);
y((y > (max(Y))) & (y < (max(Y)+buffer))) = max(Y);

z((z < (min(Z))) & (z > (min(Z)-buffer))) = min(Z);
z((z > (max(Z))) & (z < (max(Z)+buffer))) = max(Z);

%% Interpolate
for k = 1:Ncols
    Smesh.data(:, k) = interp3(Y, X, Z, squeeze(volume(:, :, :, k)),...
        y, x, z, method, extrapval);
end

%
