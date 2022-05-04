%% gethem
% Calculates mean across set of measurements.
%
%% Description
% |hem = gethem(data, info)| takes a light-level array |data| of the format
% MEAS x TIME, and using the scan metadata in |info.pairs| averages the
% shallow measurements for each wavelength present. The result is commonly
% referred to as the |hem| of a measurement set. If there is a good
% measurements logical vector present in |info.MEAS.GI|, it will be applied
% to the data; otherwise, |info.MEAS.GI| will be set to true for all
% measurements (i.e., all measurements are assumed to be good). The
% variable |hem| is output in the format WL x TIME.
% 
% |hem = gethem(data, info, sel_type, value)| allows the user to set the
% criteria for determining shallow measurements. |sel_type| can be |'r2d'|,
% |'r3d'|, or |'NN'|, corresponding to the columns of the |info.pairs|
% table, and |value| can either take the form of a two-element |[min, max]|
% vector (for |'r2d'| and |'r3d'|), or a scalar or vector containing all
% nearest neighbor numbers to be averaged. By default, this function
% averages the first nearest neighbor.
%
%% See Also
% <regcorr_help.html regcorr> | <detrend_tts_help.html detrend_tts>