%% Make_NativeSpace_4dfp
% Calculates native space of incomplete 4dfp header.
% 
%% Description
% |header_out = Make_NativeSpace_4dfp(header_in)| checks whether
% |header_in| contains the |mmppix| and |center| fields (these are not
% always present in 4dfp files). If either is absent, a default called the
% "native space" is calculated from the other fields of the volume.
% 
%% See Also
% <LoadVolumetricData_help.html LoadVolumetricData> |
% <SaveVolumetricData_help.html SaveVolumetricData>
