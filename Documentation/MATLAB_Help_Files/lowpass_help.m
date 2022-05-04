%% lowpass
% Applies a zero-phase digital lowpass filter.
%
%% Description
% |data_out = lowpass(data_in, omegaHz, frate)| takes a light-level array
% |data_in| in the MEAS x TIME format and applies to it a forward-backward
% zero-phase digital lowpass filter at a Nyquist cutoff frequency of
% |omegaHz * (2 * frate)|, returning it as |data_out|.
% 
% This function also removes the linear component of the input data.
% 
%% See Also
% <highpass_help.html highpass> | <logmean_help.html logmean> |
% <matlab:doc('filtfilt') filtfilt>