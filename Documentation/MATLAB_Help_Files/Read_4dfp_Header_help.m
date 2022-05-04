%% Read_4dfp_Header
% Reads the |.ifh| header of a |4dfp| file.
%
%% Description
% |header = Read_4dfp_Header(filename)| reads an |.ifh| text file specified by
% |filename|, containing a number of key-value pairs. The specific pairs
% are parsed and stored as fields of the output structure |header|.
%
%% See Also
% <LoadVolumetricData_help.html LoadVolumetricData>