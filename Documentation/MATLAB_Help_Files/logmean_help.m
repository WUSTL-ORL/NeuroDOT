%% logmean
% Takes the log-ratio of raw intensity data.
%
%% Description
% |data_out = logmean(data_in)| takes a light-level data array |data_in| of
% the format MEAS x TIME, and takes the negative log of each element of a
% row divided by that row's average. The result is output into |data_out|
% in the same MEAS x TIME format.
%
% The formal equation for the |logmean| operation is:
% 
%   y = -ln(x / <x>)
% 
% The MATLAB pseudocode form is:
% 
%   y(each row) = -log(x(each row) ./ mean(x(each row)));
%
% Example: If
data = [1, 10, 100; exp(1), 10*exp(1), 100*exp(1)];
%%
% then
logmean(data)
%%
% yields 
% 
%   ans = 
%           3.6109  1.3083  -.9943  3.6109  1.3083  -.9943
%
%% See Also
% <lowpass_help.html lowpass> | <highpass_help.html highpass>