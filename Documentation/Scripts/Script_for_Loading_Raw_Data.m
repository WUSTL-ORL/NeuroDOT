%% Script for Loading Raw Data
% In this script, you will be able to load raw fNIRS/DOT data to be used
% with NeuroDOT. NeuroDOT is compatible with both the SNIRF and NIRS data
% formats, in addition to NeuroDOT's native *.mat data format.
%
% Regardless of file type, after loading your raw data with NeuroDOT, you
% will have two variables: data and info.
%
% 'data' contains the raw data in the format # of channels by # of samples
% 'info' contains important information about the data that is necessary
% for processing the data with NeuroDOT.
%
% Prior to running this script, change your current directory to the 
% location of your raw data. 


%% 1. Select File to Load
% Your raw data filename including the extension (.mat, .snirf, .nirs)
filename = 'NeuroDOT_Data_Sample_CCW1.mat'; 
[pn, fn, ext] = fileparts(filename);


%% 2. Load Data
switch ext
    % Load Data in NeuroDOT (.mat) Format
    case '.mat'
        load(filename);
    
    % Load Data in SNIRF (.snirf) Format
    case '.snirf'
        [data, info] = snirf2ndot(filename, 1, [fn, '_neurodot']);
        
    % Load Data in NIRS (.nirs) Format
    case '.nirs'
        [data, info] = nirs2ndot(filename, 1, [fn, '_neurodot']);
end

