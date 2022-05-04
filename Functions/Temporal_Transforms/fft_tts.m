function [ftmag, ftdomain, ftpower, ftphase] = fft_tts(data, framerate)

% FFT_TTS Performs an FFT.
%
%   ftmag = FFT_TTS(data, framerate) takes a data array "data" of
%   the MEAS x TIME format, pads it timewise to the next highest power of
%   two (for better performance), performs the fast Fourier transform of
%   each channel using the built-in MATLAB function FFT, normalizes by the
%   padded time length, and takes the first half of the transformed data
%   (which is the positive half of the frequency domain). The result is
%   output into "ftmag".
%
%   [ftmag, ftdomain] = FFT_TTS(data, framerate) also returns the
%   corresponding normalized frequency domain "ftdomain", which extends
%   from 0 to the Nyquist frequency, which is calculated from the input
%   "framerate".
%
%   [ftmag, ftdomain, ftpower] = FFT_TTS(data, framerate) also returns the
%   power spectrum, which is the absolute value of the magnitude, squared.
%
%   [ftmag, ftdomain, ftpower, ftphase] = FFT_TTS(data, framerate) also
%   returns the phase at each frequency, as calculated by the MATLAB
%   ANGLE function.
%
% Dependencies: NORMALIZE2RANGE_TTS.
% 
% See Also: LOGMEAN, FFT, POW2, NEXTPOW2, ANGLE.
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
dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);
Ndft = pow2(nextpow2(Nt)); % Zero pack to a power of 2.
Nf = 1 + Ndft / 2;


%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Remove mean.
normdata = bsxfun(@minus,data,mean(data,2));

%% Prep data for FFT.
ftdomain = (framerate / 2) * linspace(0, 1, Nf); % domain of FFT: [zero:Nyquist]

%% Perform FFT.
P = fft(normdata, Ndft, 2) / Ndft; % Do FFT in TIME dimension and normalize by Ndft.
ftmag = sqrt(2) * abs(P(:, 1:Nf)); % Take positive frequencies, x2 for negative frequencies.

%% N-D Output.
if NDtf
    ftmag = reshape(ftmag, [dims(1:end-1), Nf]);
end

%% Other outputs.
ftpower = abs(ftmag) .^ 2;
ftphase = angle(P(:, 1:Nf));



%
