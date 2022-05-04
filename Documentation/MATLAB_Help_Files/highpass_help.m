%% highpass
% Applies a zero-phase digital highpass filter.
%
%% Description
% |data_out = highpass(data_in, omegaHz, frate)| takes a light-level array
% |data_in| in the MEAS x TIME format and applies to it a forward-backward
% zero-phase digital highpass filter at a Nyquist cutoff frequency of
% |omegaHz * (2 * frate)|, returning it as |data_out|.
% 
% This function also removes the linear component of the input data.
% 
%% See Also
% <lowpass_help.html lowpass> | <logmean_help.html logmean> |
% <matlab:doc('filtfilt') filtfilt>
