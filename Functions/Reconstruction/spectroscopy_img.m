function cortex_hb = spectroscopy_img(cortex_mu_a, E)

% SPECTROSCOPY_IMG Completes the Beer-Lambert law from a reconstructed
% image.
%
%   img_out = SPECTROSCOPY_IMG(img_in, E) takes the reconstructed VOX x
%   TIME x WL mua image "img_in" and multiplies it by the inverse
%   of the extinction coefficient matrix "E" to create an output VOX x TIME
%   x HB matrix "img_out", where HB 1 and 2 are the voxel-space time series
%   images for HbO and HbR, respectively.
%
% See Also: RECONSTRUCT_IMG, AFFINE3D_IMG.
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
[Nvox, Nt, Nc] = size(cortex_mu_a);
[E1, E2] = size(E);
umol_scale = 1000;

%% Check compatibility of the image and E matrix.
if E1 ~= Nc
    error(['ERROR: The image wavelengths and spectroscopy matrix dimensions do not match.'])
end

%% Invert Spectroscopy Matrix.
iE = inv(E);

% Initialize Outputs
cortex_hb = zeros(Nvox, Nt, Nc, 'single');
for k = 1:E2
    temp = zeros(Nvox, Nt);
    for l = 1:E1
        temp = temp + squeeze(iE(k, l)) .* squeeze(cortex_mu_a(:, :, l));
    end
    cortex_hb(:, :, k) = temp;
end

cortex_hb = cortex_hb .* umol_scale; % Fix units to umol



%
