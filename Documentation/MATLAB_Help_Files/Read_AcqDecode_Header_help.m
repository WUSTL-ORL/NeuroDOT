%% Read_AcqDecode_Header
% Reads the header section of a binary NIRS file.
% 
%% Description
% |[Nd, Ns, Nwl, Nt] = Read_AcqDecode_Header(fid)| reads the first bytes of
% an open binary NIRS file identified by |fid|, and selects the correct
% format to load the rest of the header.
% 
%% See Also
% <LoadMulti_AcqDecode_Data_help.html LoadMulti_AcqDecode_Data> |
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data> |
% <matlab:doc('fopen') fopen> | <matlab:doc('fread') fread>