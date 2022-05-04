%% resample_tts
% Resample data while maintaining linear signal component.
%
%% Description
% |[data_out, info_out] = resample_tts(data_in, info_in, omega_resample, tol,
% framerate)| takes a raw light-level data array |data_in| of the format
% MEAS x TIME, and resamples it (typically downward) to a new frequency
% using the built-in MATLAB function |resample|. The new sampling frequency
% is calculated as the ratio of input |omega_resample| divided by
% |framerate| (both scalars), to within the tolerance specified by |tol|.
%
% This function is needed because the linear signal components, which can
% be important in other NeuroDOT pipeline calculations, can be
% inadvertently removed by downsampling using |resample| alone.
%
% Note: This function resamples synch points in addition to data. Be sure
% to take care that your data and synch points match after running this
% function! |info.paradigm.init_synchpts| stores the original synch points
% if you need to restore them.
% 
%% See Also
% <detrend_tts_help.html detrend_tts> | <matlab:doc('resample') resample>