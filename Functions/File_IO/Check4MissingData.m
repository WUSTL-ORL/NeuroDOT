function [data_out, info_out, framepts_out] = Check4MissingData(data_in, info_in, framepts_in)

% CHECK4MISSINGDATA Validates framesynch file against key file info.
%
%   [data_out, info_out, framepts_out] = CHECK4MISSINGDATA(data_in,
%   info_in, framepts_in) takes the -info.txt key file information and
%   compares it to the .fs framesynch file. If any frames are missing
%   samples, the user is asked to abort, exclude incomplete frame(s), or
%   cut all samples after the first incomplete frame.
%
% See Also: LOAD_ACQDECODE_DATA.
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

%% Parameters and initialization.
data_out = data_in;
info_out = info_in;
framepts_out = framepts_in;

% Approximate frame length = (time steps) * (flashing time + buffer time)
fguess = info_out.io.nts * (info_out.io.nsamp + info_out.io.nblank);

% Actual frame lengths.
fsizes = diff(framepts_out);

%% Check for missing data
bad_frames = (fsizes < fguess - 3)  |  (fsizes > fguess + 3);

%% Discuss with user.
if any(bad_frames)
    figure
    plot(fsizes, 'ok')
    ylabel('Frame Size')
    xlabel('Frame Number')
    disp('** The frame'' sizes exceed expected deviation: samples likely missing **')
    switch input('1: Abort, 2: Exclude All, 3: Keep Previous')
        case 1
            error('** Aborted due to error in data recording **')
        case 2
            data_out = data_in(:, :, :, ~bad_frames); % Note: This is significantly different from the ND1 version, but does the same thing.
            framepts_out = framepts_out(~bad_frames);
            disp(['The following frames were excluded: ', num2str(find(bad_frames))])
        case 3
            stoppt = find(bad_frames, 'first') - 2;
            data_out = data_in(:, :, :, 1:stoppt);
            framepts_out = framepts_out(1:stoppt);
            fsizes = fsizes(1:stoppt);
            disp(['The data was terminated after frame: ', num2str(stoppt)])
    end
end

%% Save information.
info_out.io.framesize = mean(fsizes); % mean frame length



%
