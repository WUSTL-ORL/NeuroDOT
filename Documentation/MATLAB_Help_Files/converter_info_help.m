%% converter_info 
% Converts the ND1 info structure to ND2 format.
% 
%% Description
% |info_out = converter_info(info_in, conversion)| adapts the metadata
% structures between ND1 and ND2.
% 
% |conversion| can be either |'ND2 to ND1'| or |'ND1 to ND2'|.
% 
% All leftover fields from |'ND1 to ND2'| are stored in |info.misc|.
% 
% Only intended for use with ND1 and ND2 data structures.
% 
%% See Also
% <converter_data_help.html converter_data> |
% <converter_ND1_to_ND2_help.html converter_ND1_to_ND2> |
% <converter_ND2_to_ND1_help.html converter_ND2_to_ND1>