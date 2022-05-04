function iA_out = smooth_Amat(iA_in, dim, gsigma)

% SMOOTH_AMAT Performs Gaussian smoothing on a sensitivity matrix.
%
%   iA_out = SMOOTH_AMAT(iA_in, dim, gsigma) takes the inverted VOX x MEAS
%   sensitivity matrix "iA_in" and performs Gaussian smoothing in the 3D
%   voxel space on each of the concatenated measurement matrices within,
%   returning it as "iA_out". The user specifies the Gaussian filter
%   half-width "gsigma".
%
% See Also: TIKHONOV_INVERT_AMAT, RECONSTRUCT_IMG, FINDGOODMEAS.
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
[Nvox, Nm] = size(iA_in);
iA_out = zeros(Nvox, Nm, 'single');

nVx = dim.nVx;
nVy = dim.nVy;
nVz = dim.nVz;

gsigma = gsigma / dim.sV;


%% Preallocate voxel space.
if isfield(dim, 'Good_Vox')
    GV = dim.Good_Vox;
else
    GV = ones(nVx, nVy, nVz); % WARNING: THIS RUNS WAY SLOWER.
end

%% Do smoothing in parallel.
parfor k = 1:Nm
    iAvox = zeros(nVx, nVy, nVz); % Set up temp iAvox
    iAvox(GV) = iA_in(:, k); % Grab iA vox for a meas
    iAvox = imgaussfilt3(iAvox,gsigma); % smooth
    iA_out(:, k) = single(iAvox(GV)); % Put back into vector form on the way out.
end


%
