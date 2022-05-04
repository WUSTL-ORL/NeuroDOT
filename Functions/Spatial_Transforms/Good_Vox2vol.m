function imgvol = Good_Vox2vol(img, dim)

% GOOD_VOX2VOL Turns a Good Voxels data stream into a volume.
% 
%   imgvol = GOOD_VOX2VOL(img, dim) reshapes a VOX x TIME array "img" into
%   an X x Y x Z x TIME array "imgvol", according to the dimensions of the
%   space described by "dim".
% 
% See Also: SPECTROSCOPY_IMG.
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
[Nvox, Nt] = size(img);
if Nvox==1 && Nt>1
    img=img';
    [Nvox, Nt] = size(img);
end

%% Stream image into good voxels.
imgvol = zeros(dim.nVx * dim.nVy * dim.nVz, Nt);
imgvol(dim.Good_Vox, :) = img;

%% Reshape image into the voxel space.
imgvol = reshape(imgvol, dim.nVx, dim.nVy, dim.nVz, Nt);