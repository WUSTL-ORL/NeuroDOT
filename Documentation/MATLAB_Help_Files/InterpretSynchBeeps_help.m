%% InterpretSynchBeeps
% Automated synch file analysis.
% 
%% Description
% |[synchpts, synchtype] = InterpretSynchBeeps(synch)| uses the MATLAB
% function |findpeaks| to convert the stim synch array |synch| into a set
% of pulse timepoints and pulse heights, which are returned in |synchpts|
% and |synchtype| respectively.
% 
%% See Also
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data> | <ReadAux_help.html
% ReadAux> | <InterpretStimSynch_help.html InterpretStimSynch> |
% <matlab:doc('findpeaks') findpeaks>