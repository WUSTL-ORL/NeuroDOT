%% converter_ND2_to_HOMER
% Converts workspace from ND2 to HOMER.
%
%% Description
% |nirs = converter_ND2_to_HOMER(data, info, aux)| takes MEAS x TIME array
% |data|, structure |info|, and optionally the |aux| structure, and
% converts them to a |nirs| structure that can be saved into a MAT file
% with the |.nirs| extension, and is compatible with HOMER2. 
%  
%% See Also
% <converter_HOMER_to_ND2_help.html converter_HOMER_to_ND2>