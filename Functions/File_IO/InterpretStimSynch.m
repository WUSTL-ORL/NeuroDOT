function synch = InterpretStimSynch(raw_synch, framepts, info)

% INTERPRETSTIMSYNCH Reads stimulus from frame synch file.
% 
%    synch = INTERPRETSTIMSYNCH(raw_synch, framepts, info) reads the frames
%    "framepts" of an auxiliary synch file "raw_synch", then calculates and
%    returns synch block standard deviation and frequency information in
%    "synch".
% 
% See Also: LOAD_ACQDECODE_DATA, READAUX, INTERPRETSYNCHBEEPS.
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
synch = zeros(2, info.io.nframe);

%% Cycle over all frames in recording
for f = 1:info.io.nframe - 1
    % Pick out a frame from synch channel.
    synchblock = raw_synch(framepts(f):(framepts(f) + floor(info.io.framesize) - 1));
    
    % Take std to look for beeps.
    synch(1, f) = std(synchblock);
    
    % Find frequency at peak.
    fsynch = abs(fft(synchblock));
    [~, synch(2, f)] = max(fsynch(1:floor(info.io.framesize / 2)));
end



%
