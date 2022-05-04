%% normalize2range_tts
% Normalizes time traces to average of a range.
% 
%% Description
% |data_out = normalize2range_tts(data_in, range_to_avg)| takes a light-level
% array |data_in| in the MEAS x TIME format, and calculates the mean of
% each channel for the time points specified by vector |range_to_avg|. If
% no range is specified, the function defaults to averaging all time
% points. Each channel then has this mean subtracted from it, and the
% result is output as |data_out|.
% 
% The formal equation for the |normalize2range_tts| operation is:
% 
%    y = x - <x(range)>
% 
% The MATLAB pseudocode form is:
% 
%    y(each row) = x(each row) - mean(x(range));
%
% Example 1: If
data = [21, 20, 19; 300, 295, 305];
%%
% then
normalize2range_tts(data)
%% 
% yields
% 
%   ans = 
%           1   0   -1
%           0   -5  5
%  
% Example 2: If
data = [1, 2, 3, 4, 5, 6, 7];
%%
% then
normalize2range_tts(data, 1:5)
%%
% yields
% 
%   ans = 
%           -2  -1  0   1   2   3   4
% 
%% See Also
% <fft_tts_help.html fft_tts> | <logmean_help.html logmean>