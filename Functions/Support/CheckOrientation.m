function [xdir, ydir, zdir] = CheckOrientation(info)

% CHECKORIENTATION Interpretes 3D orientation information.
%
%   [xdir, ydir, zdir] = CHECKORIENTATION(info) checks the directions of
%   the data encoding in "info.optodes.plot3orientation" to make it
%   compatible with 3D cap visualizations. The default reverses the x-axis
%   and leaves the other two alone. The directions are returned in "xdir",
%   "ydir", and "zdir" 
%
% THIS FUNCTION IS BEING DEPRECATED. Until it is cleaned out, its
% functionality is being steralized.
% 
% See Also: PLOTCAPDATA.
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
xdir = 'normal';
ydir = 'normal';
zdir = 'normal';

% %% Interpret the orientation codes.
% switch info.optodes.plot3orientation.i
%     case 'R2L'
%         xdir = 'reverse';
%     case 'L2R'
%         xdir = 'normal';
%     otherwise
%         xdir = 'reverse'; % Default.
% end
% switch info.optodes.plot3orientation.j
%     case 'A2P'
%         ydir = 'reverse';
%     case 'P2A'
%         ydir = 'normal';
%     otherwise
%         ydir = 'normal'; % Default.
% end
% switch info.optodes.plot3orientation.k
%     case {'V2D', 'I2S'}
%         zdir = 'reverse';
%     case {'D2V', 'S2I'}
%         zdir = 'normal';
%     otherwise
%         zdir = 'normal'; % Default.
% end