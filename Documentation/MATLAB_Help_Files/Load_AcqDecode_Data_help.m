%% Load_AcqDecode_Data
% Loads a single AcqDecode file.
%
%% Description
% |[data, info, synch, aux, framepts] = Load_AcqDecode_Data(filename)|
% reads a raw data file and all secondary files scan files in the directory
% specified by |filename|. The raw data is validated and checked for
% missing frames, then reshaped into the ND2 format.
%
% All secondary files must be located in the same folder as |filename|.
%
% The information is returned as a raw light level array |data|, a metadata
% structure |info|, stimulus synchronization information |synch|, raw
% auxiliary file data |aux|, and frame time point array |framepts|.
%
% Currently supports |.mag| and |.iq| files.
% 
%% Dependencies
% <Read_AcqDecode_Header_help.html Read_AcqDecode_Header> |
% <ReadAux_help.html ReadAux> | <ReadInfoTxt_help.html ReadInfoTxt> |
% <Check4MissingData_help.html Check4MissingData> |
% <InterpretStimSynch_help.html InterpretStimSynch> |
% <InterpretSynchBeeps_help.html InterpretSynchBeeps> |
% <InterpretPulses_help.html InterpretPulses>
% 
%% See Also
% <LoadMulti_AcqDecode_Data_help.html LoadMulti_AcqDecode_Data>