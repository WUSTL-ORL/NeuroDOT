function [IM2, info2] = Zero_Pad_Vol(IM,info1,Np)
%
% This function takes an image volume and adds zeros at each edge.
% For now, the volume is assumed to be 3D and the padding is uniform for
% all 6 sides (padding = Np voxels).

% Inputs
% IM: image volume
% info1: metadata for image volume
% Np: number of pixels to pad image volume by on each side

% Outputs
% IM2: padded image volume
% info2: metadata for padded image volume

% Copyright (c) 2017 Washington University 
% Created By: Ari Segel
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

%% Copy info
info2 = info1;

%% Edit info structure with padding info
% Dimensions first
info2.nVx = info2.nVx + 2*Np;
info2.nVy = info2.nVy + 2*Np;
info2.nVz = info2.nVz + 2*Np;

% Re-calculate center
info2.center = info1.center + Np .* info1.mmppix;

IM2 = zeros(size(IM) + 2*Np);
IM2((Np+1):(end-Np), (Np+1):(end-Np),...
    (Np+1):(end-Np)) = IM;


end 
