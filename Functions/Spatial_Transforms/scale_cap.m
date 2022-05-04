function tpos_out = scale_cap(tpos_in, M)

% SCALE_CAP Scales a cap grid.
% 
%   tpos_out = SCALE_CAP(tpos_in, M) scales the full cap grid given in
%   "tpos_in" by a factor of "M" around its centroid, and outputs it as
%   "tpos_out". If M is a 3-component vector, then the cap is scaled about
%   its centroid by a factor of M(1) in x, M(2) in y, and M(3) in z.
% 
% See Also: CAP_FITTER.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht, Zachary E. Markow
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
if numel(M) == 1
    M_new = repmat(M,3,1);
else
    M_new = M;
end

centroid = mean(tpos_in, 1);
centroid_mat = repmat(centroid, [size(tpos_in, 1), 1]);

scaling_mat = diag(M_new);

%% Do Scaling.
tpos_out = (tpos_in - centroid_mat) * scaling_mat + centroid_mat;



%
