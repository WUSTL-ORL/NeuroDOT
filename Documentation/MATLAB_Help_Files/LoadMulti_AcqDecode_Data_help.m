%% LoadMulti_AcqDecode_Data
% Loads and combines data from multiple AcqDecode files for a single scan.
%
%% Description
% |[data, info, synch, aux] = LoadMulti_AcqDecode_Data(filename, path)
% finds any acquisitions |filename| located in directory |path|, and loads
% the raw data files and any other relevant scan files contained therein.
% The information is returned as a raw light level array |data|, a metadata
% structure |info|, stimulus synchronization information |synch|, and raw
% auxiliary file data |aux|.
% 
% It is assumed that all files use the AcqDecode format:
% |"DATE-SUBJECT-TAGletter"|, with a trailing lowercase letter for multiple
% acquisitions. Also, separate acquisitions are stored in subfolders named
% by |"DATEletter"|, EG |"150115a"| and |"150115b"|.
%
% Currently supports |.mag| and |.iq| files.
% 
% |[data, info, synch, aux] = LoadMulti_AcqDecode_Data(filename, path,
% flags)| uses the |flags| structure to specify loading parameters.
% 
%% Parameters
% |flags| fields that apply to this function (and their defaults):
% 
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>Nsys</td><td>2</td><td>Number of acquisitions.
% </table></html>
% 
%% Dependencies
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data> |
% <Crop2Synch_help.html Crop2Synch>