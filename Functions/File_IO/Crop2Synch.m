function [data_out, info_out] = Crop2Synch(data_in, info_in, flags)

% CROP2SYNCH Crops data to synch points.
% 
%   [data_out, info_out] = CROP2SYNCH(data_in, info_in, flags) crops a MEAS
%   x TIME array "data_in" to the pulses in "info.paradigm.synchpts",
%   adjusts the synch points to the new scale, saves the original synch
%   points, and returns both structures. This function uses the "flags"
%   structure to specify loading parameters.
% 
%   "flags" fields that apply to this function (and their defaults):
%       crop_level      (0 or 2)    Data cropping. 0 = none, 1 = start
%                                   pulse, 2 = start and stop.
% 
% See Also: LOAD_ACQDECODE_DATA, INTERPRETSYNCHBEEPS.
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
data_out = data_in;
info_out = info_in;

if ~exist('flags', 'var')
    flags = [];
end

if ~isfield(flags, 'crop_level')
    if isfield(flags, 'Nsys')  &&  flags.Nsys == 1
        flags.crop_level = 0;
    else
        flags.crop_level = 2;
    end
end

%% Do cropping.
if numel(info_out.paradigm.synchpts) > 0
    if flags.crop_level >= 1
        % Crop to start pulse.
        data_out(:, 1:info_out.paradigm.synchpts(1)-1) = [];
        
        % Adjust synch points.
        info_out.paradigm.synchpts = info_out.paradigm.synchpts - info_out.paradigm.synchpts(1) + 1;
        
        % Save these as the initial synch pts.
        info_out.paradigm.init_synchpts = info_out.paradigm.synchpts;
    end
    if flags.crop_level == 2
        data_out(:, info_out.paradigm.synchpts(end)+1:end) = []; % Crop to stop pulse.
    end
end



%
