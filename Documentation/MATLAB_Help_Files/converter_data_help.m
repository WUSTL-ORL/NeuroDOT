%% converter_data
% Reshapes data between ND1 and ND2 formats.
% 
%% Description
% |data_out = converter_data(data_in)| reshapes any ND1 array |data_in| to
% ND2 array format in |data_out|, or vice versa. The direction is
% determined by the string |conversion|, and in the |'ND2 to ND1'|
% direction, |info| is needed.
% 
% |conversion| can be either |'ND2 to ND1'| or |'ND1 to ND2'|.
% 
%% See Also
% <converter_ND1_to_ND2_help.html converter_ND1_to_ND2> |
% <converter_ND2_to_ND1_help.html converter_ND2_to_ND1> |
% <converter_info_help.html converter_info>