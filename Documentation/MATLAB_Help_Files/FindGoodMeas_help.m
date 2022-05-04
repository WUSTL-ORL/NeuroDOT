%% FindGoodMeas
% Performs "Good Measurements" analysis.
%
%% Description
% |info_out = FindGoodMeas(data, info_in)| takes a light-level array |data|
% in the MEAS x TIME format, and calculates the std of each channel as
% its noise level. These are then thresholded by the default value of
% |0.075| to create a logical array, and both the are returned as MEAS x 1
% columns of the |info.MEAS| table. If pulse synch point information exists
% in |info.system.synchpts|, then FindGoodMeas will crop the data to the
% start and stop pulses.
% 
% |info_out = FindGoodMeas(data, info_in, bthresh)| allows the user to
% specify a threshold value.
%
%% See Also
% <PlotCapGoodMeas_help.html PlotCapGoodMeas> | <PlotHistogramSTD_help.html
% PlotHistogramSTD>