%% fft_tts
% Performs an FFT.
%
%% Description
% |ftmag = fft_tts(data, framerate)| takes a data array |data| of
% the MEAS x TIME format, pads it timewise to the next highest power of
% two (for better performance), performs the fast Fourier transform of
% each channel using the built-in MATLAB function |fft|, normalizes by the
% padded time length, and takes the first half of the transformed data
% (which is the positive half of the frequency domain). The result is
% output into |ftmag|.
% 
% |[ftmag, ftdomain] = fft_tts(data, framerate)| also returns the
% corresponding normalized frequency domain |ftdomain|, which extends
% from 0 to the Nyquist frequency, which is calculated from the input
% |framerate|.
% 
% |[ftmag, ftdomain, ftpower] = fft_tts(data, framerate)| also returns the
% power spectrum, which is the absolute value of the magnitude, squared.
% 
% |[ftmag, ftdomain, ftpower, ftphase] = fft_tts(data, framerate)| also
% returns the phase at each frequency, as calculated by the MATLAB
% |angle| function.
% 
%% Dependencies
% <normalize2range_tts_help.html normalize2range_tts>
% 
%% See Also
% <logmean_help.html logmean> | <matlab:doc('fft') fft> |
% <matlab:doc('pow2') pow2> | <matlab:doc('nextpow2') nextpow2> |
% <matlab:doc('angle') angle>