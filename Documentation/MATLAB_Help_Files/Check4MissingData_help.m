%% Check4MissingData
% Validates framesynch file against key file info.
%
%% Description
% |[data_out, info_out, framepts_out] = Check4MissingData(data_in, info_in,
% framepts_in)| takes the |-info.txt| key file information and compares it
% to the |.fs| framesynch file. If any frames are missing samples, the user
% is asked to abort, exclude incomplete frame(s), or cut all samples after
% the first incomplete frame.
% 
%% See Also
% <Load_AcqDecode_Data_help.html Load_AcqDecode_Data>