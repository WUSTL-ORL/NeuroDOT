%% ReadInfoTxt
% Reads a key file.
% 
%% Description
% |io = ReadInfoTxt(filename)| scans the |-info.txt| key file specified by
% |filename| and parses its contents into key-value pairs, which are
% returned in the |io| substructure of the main ND2 |info| structure. Most
% of the pairs are either only used in the loading process or kept for
% later troubleshooting if the data set shows inconsistencies.
% 
%% See Also
% <LoadMulti_AcqDecode_Data_help.html LoadMulti_AcqDecode_Data> |
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data> |
% <ReadAux_help.html ReadAux>