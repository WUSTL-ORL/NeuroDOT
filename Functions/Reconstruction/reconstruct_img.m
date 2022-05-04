function cortex_mu_a = reconstruct_img(lmdata, iA)

% RECONSTRUCT_IMG Performs image reconstruction by wavelength using
% inverted A-matrix.
%
%   img = RECONSTRUCT_IMG(data, iA) takes the inverted VOX x MEAS
%   sensitivity matrix "iA" and right-multiplies it by the logmeaned MEAS x
%   TIME light-level matrix "lmdata" to reconstruct an image in voxel
%   space. The image is output in a VOX x TIME matrix "img".
%
% See Also: TIKHONOV_INVERT_AMAT, SMOOTH_AMAT, SPECTROSCOPY_IMG,
% FINDGOODMEAS.
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
units_scaling = 1/100; % Assuming OptProp in mm^-1 %% Will move to 10 shortly.

%% Reconstruct.
cortex_mu_a = iA * lmdata;

%% Correct units and convert to single.
cortex_mu_a = single(cortex_mu_a .* units_scaling);