%% Crop2Synch
% Crops data to synch points.
% 
%% Description
% |[data_out, info_out] = Crop2Synch(data_in, info_in, flags)| crops a MEAS
% x TIME array |data_in| to the pulses in |info.paradigm.synchpts|, adjusts
% the synch points to the new scale, saves the original synch points, and
% returns both structures. This function uses the |flags| structure to
% specify loading parameters.
% 
%% Parameters
% |flags| fields that apply to this function (and their defaults):
% 
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>crop_level</td><td>(0 or 2)</td><td>Data cropping. 0 = none, 1 =
% start pulse, 2 = start and stop.</td></tr>
% </table></html>
% 
%% See Also
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data> |
% <InterpretSynchBeeps_help.html InterpretSynchBeeps>