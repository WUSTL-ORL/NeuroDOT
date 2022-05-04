function header_out = Make_NativeSpace_4dfp(header_in)

% MAKE_NATIVESPACE_4DFP Calculates native space of incomplete 4dfp header.
% 
%   header_out = MAKE_NATIVESPACE_4DFP(header_in) checks whether
%   "header_in" contains the "mmppix" and "center" fields (these are not
%   always present in 4dfp files). If either is absent, a default called
%   the "native space" is calculated from the other fields of the volume.
% 
% See Also: LOADVOLUMETRICDATA, SAVEVOLUMETRICDATA.
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
header_out = header_in;

%% Check for mmppix, calculate if not present.
if ~isfield(header_out, 'mmppix')
    if all(~isfield(header_out, {'mmx', 'mmy', 'mmz'}))
        error(['*** Error: Required field(s) ''mmx'', ''mmy'', or ''mmz''',...
            ' not present in input header. ***'])
    end
    
    header_out.mmppix = [header_out.mmx, -header_out.mmy, -header_out.mmz];
    
end

%% Check for center, calculate if not present.
if ~isfield(header_out, 'center')
    if all(~isfield(header_out, {'nVx', 'nVy', 'nVz'}))
        error(['*** Error: Required field(s) ''nVx'', ''nVy'', or ''nVz''',...
            ' not present in input header. ***'])
    end
    
    header_out.center = header_out.mmppix .* ( [header_out.nVx,...
        header_out.nVy, header_out.nVz] / 2 + [0, 1, 1]);
    
end



%
