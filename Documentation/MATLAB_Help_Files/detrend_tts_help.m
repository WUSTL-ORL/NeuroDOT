%% detrend_tts
% Performs linear detrending.
%
%% Description 
% |data_out = detrend_tts(data_in)| takes a raw light-level data array
% |data_in| of the format MEAS x TIME and removes the straight-line fit
% along the TIME dimension from each measurement, returning it as
% |data_out|.
% 
%% See Also
% <logmean_help.html logmean>