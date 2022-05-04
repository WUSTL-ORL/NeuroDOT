function [synchpts, synchtype] = InterpretSynchBeeps(synch)

% INTERPRETSYNCHBEEPS Automated synch file analysis.
% 
%   [synchpts, synchtype] = INTERPRETSYNCHBEEPS(synch) uses the MATLAB
%   function FINDPEAKS to convert the stim synch array "synch" into a set
%   of pulse timepoints and pulse heights, which are returned in "synchpts"
%   and "synchtype" respectively.
% 
% See Also: LOAD_ACQDECODE_DATA, READAUX, INTERPRETSTIMSYNCH, FINDPEAKS.
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
% The 1st row contains the std of signal intensity in each frame.
synch1 = squeeze(synch(1, :)) - mean(synch(1, :)); 

% Set peak threshold.
thresh = 0.25;
m = thresh * max(synch1);
dt = 5;

jitter_factor = 5;

%% Peak finding.
[~, synchpts] = findpeaks(synch1, 'minpeakheight', m, 'minpeakdistance', dt);

% If we have frequency info (2nd row), then look at peak times.
if size(synch, 1) > 1
    synchtype = synch(2, synchpts);
    
    % Fix any frequency space jitter by quantizing to 5.
    synchtype = jitter_factor * round(synchtype / jitter_factor);
else
    synchtype = [];
end



%
