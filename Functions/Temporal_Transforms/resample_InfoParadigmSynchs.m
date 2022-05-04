function [info_out] = resample_InfoParadigmSynchs(info_in, omega_resample, tol, framerate)

% resample_InfoParadigmSynchs Resample the timing of the synchpoints
%   within info.paradigm.synchpts to remain consistent with data.
%
%   [info_out] = RESAMPLE_INFOPARADIGMSYNCHS(info_in, tHz, tol,framerate) 
%   takes info structure containing info.paradigm that includes the timing
%   of the experiment with synchronization event times listed in
%   info.paradigm.synchpts. These are resampled to a new frequency
%   using the built-in MATLAB function RESAMPLE. The new sampling frequency
%   is calculated as the ratio of input "omega_resample" divided by
%   "framerate" (both scalars), to within the tolerance specified by "tol".
%
%
% See Also: DETREND_TTS, RESAMPLE.
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
info_out = info_in;

if ~exist('framerate', 'var')  ||  isempty(framerate)
    if isfield(info_in, 'system')  &&  ~isempty(info_in.system)...
            && isfield(info_in.system, 'framerate') ...
            && ~isempty(info_in.system.framerate)
        framerate = info_in.system.framerate;
    end
end
info_out.system.init_framerate=framerate;

if ~exist('omega_resample', 'var')  ||  isempty(omega_resample)
    omega_resample = 1;
end
if ~exist('tol', 'var')  ||  isempty(tol)
    tol = 1e-5;
end


%% Approximate desired resampling ratio as a fraction.
[N, D] = rat(omega_resample / framerate, tol);
info_out.system.framerate = omega_resample;


%% Fix synch pts to new framerate.
if isfield(info_in,'paradigm')
    if isfield(info_in.paradigm, 'init_synchpts')
        info_out.paradigm.synchpts = round(N .* info_out.paradigm.init_synchpts ./ D);
        info_out.paradigm.synchpts(info_out.paradigm.synchpts == 0) = 1;
    else        
        info_out.paradigm.init_synchpts=info_out.paradigm.synchpts;
        info_out.paradigm.synchpts = round(N .* info_out.paradigm.synchpts ./ D);
        info_out.paradigm.synchpts(info_out.paradigm.synchpts == 0) = 1;
    end
end


%
