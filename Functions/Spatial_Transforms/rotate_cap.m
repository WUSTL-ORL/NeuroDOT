function tpos_out = rotate_cap(tpos_in, dTheta)

% ROTATE_CAP Rotates the cap in space.
%
%   tpos_out = ROTATE_CAP(tpos_in, dTheta) rotates the cap grid given by
%   "tpos_in" by the rotation vector "dTheta" (in degrees) and outputs it
%   as "tpos_out".
% 
% Dependencies: ROTATION_MATRIX.
% 
% See Also: PLOTLRMESHES, SCALE_CAP.
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
centroid = mean(tpos_in, 1);
centroid_mat = repmat(centroid, [size(tpos_in, 1), 1]);
Dx_axis = dTheta(1); % Pos rotates true right down
Dy_axis = dTheta(2); % Pos rotates CCW from top
Dz_axis = dTheta(3); % Pos rotates back of cap up
d2r = pi / 180; % Convert from degrees to radians

%% Create rotation matrices.
rotX = rotation_matrix('x', Dx_axis * d2r);
rotY = rotation_matrix('y', Dy_axis * d2r);
rotZ = rotation_matrix('z', Dz_axis * d2r);

%% Rotate around Centroid.
rot = rotX * rotY * rotZ;
tpos_out = (tpos_in - centroid_mat) * rot + centroid_mat;



%
