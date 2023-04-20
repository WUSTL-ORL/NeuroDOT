%% Script for Generating a pad file

%% NOTE
% When running this script, please DO NOT HIT THE GREEN TRIANGLE "RUN" BUTTON
% Please run individual code sections by highlighting and evalutating the highlighted selection
% Or, you can hit "ctrl + enter" on each section of code once the section's background turns beige

%% Pathing and output directory
% Please make sure you have the neuroDOT toolbox on your matlab path

% Set output directory
outputdir = (''); %path to directory goes in between quotes
cd(outputdir)

%% Load data and get optode locations

% load data that contains optode locations and wavelengths
Data = load('NeuroDOT/Data/NeuroDOT_Data_Sample_CCW1');

% create grid structure and place optode locations inside
    % grid structure should contain the following fields at a minimum
    % 3D source positions: grid.spos3
    % 3D detector positions: grid.dpos3
grid = struct;
grid.spos3 = Data.info.optodes.spos3;
grid.dpos3 = Data.info.optodes.dpos3;

% If available, place 2D optode positions in grid structure 
grid.spos2 = Data.info.optodes.spos2;
grid.dpos2 = Data.info.optodes.dpos2;


%% Create info structure (this is the pad file)

params.lambda = unique(Data.info.pairs.lambda); % Ex: [750, 850]
params.mod = 'CW'; % Modulation type or frequency
params.CapName = 'ND_tutorial_pad'; % Create this yourself
info = Generate_pad_from_grid(grid,params);


%% Visualize layout of pad fiole

% 2D layout - PlotCap
PlotCap(info)

% 3D layout - PlotCap
params_cap.dimension = '3D';
PlotCap(info, params_cap);view([-40,30])

% 3D layout - Draw_Foci
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3); %SD positions
Ns=size(info.optodes.spos3,1); %Number of sources
Nd=size(info.optodes.dpos3,1); %Number of detectors
paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1)); % colors for sources (red) and detectors (blue)
paramsFoci.color(1,:) = [1 0.4 0.6]; %Pink for s1
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %Light blue for d1
figure;Draw_Foci_191203(tpos,paramsFoci);view([-40,30])


%% Make sanity check plots

% Visualize r2d by itself (hist)
figure;histogram(info.pairs.r2d,1000);xlabel('R_S_D');...
    ylabel('N_m_e_a_s');title('2D SD Separations');xlim([0 60])

% Visualize r3d by itself (hist)
figure;histogram(info.pairs.r3d,1000);xlabel('R_S_D');...
    ylabel('N_m_e_a_s');title('3D SD Separations');xlim([0 60])

% Visualize r3d by NN (scatter)
figure;scatter(info.pairs.r3d,info.pairs.NN);xlabel('R3d');ylabel('NN');
title([info.optodes.CapName, ', 3D SD sep by NN'], 'interpreter', 'none');


%% Save pad file
save(['Pad_', info.optodes.CapName, '.mat'], 'info')


close all; clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Match pad to data Section 
% Some data uses a different bookeeping strategy than NeuroDOT
% As such, it is important to make sure that whatever we make using
% NeuroDOT reflects the underlying data that we want to analyze
% This section will re-organize the pad file measurement list found in
% info.pairs to match the measurement list seen in the data

%% Load data file
Data = load('NeuroDOT/Data/NeuroDOT_Data_Sample_GV1');


%% Generate pad file

% Optode positions
grid = struct;
grid.spos3 = Data.info.optodes.spos3;
grid.dpos3 = Data.info.optodes.dpos3; 
grid.spos2 = Data.info.optodes.spos2;
grid.dpos2 = Data.info.optodes.dpos2;

% Params
params.lambda = unique(Data.info.pairs.lambda); % Ex: [750, 850]
params.mod = 'CW'; % Modulation type or frequency
params.CapName = 'ND_tutorial_pad_cropSection'; % Create this yourself

% Generate
info = Generate_pad_from_grid(grid,params);

% Visualize 2D
PlotCap(info)


%% Get Measurement lists from Pad and from Data

% Create temp pairs structure, this will get modified to match data
temp = info.pairs;

% Make measurement list from pad
pad_measList_before = [temp.Src, temp.Det, temp.WL];

% Make measurement list from data
data_measList = [Data.info.pairs.Src, Data.info.pairs.Det,...
    Data.info.pairs.WL];

% Are measurement lists in the same order?
order_before = isequal(data_measList, pad_measList_before)


%% Crop to match data
% Get order of measmurements in data
[~,Idx]=ismember(data_measList,pad_measList_before,'rows');

% Re-order pairs structure
temp.r3d=[temp.r3d(Idx)];
temp.r2d=[temp.r2d(Idx)];
% temp.r2d = temp.r3d; %for sparse pads, uncomment this line - set r2d = to r3d
temp.NN = [temp.NN(Idx)];
temp.Src = [temp.Src(Idx)];
temp.Det = [temp.Det(Idx)];
temp.WL = [temp.WL(Idx);]; %Get correct size
temp.lambda = [temp.lambda(Idx)]; %Get corect size
temp.Mod = [temp.Mod(Idx,:),];

% Make sure measurements are identical
pad_measList_after = [temp.Src, temp.Det, temp.WL];
order_after = isequal(data_measList, pad_measList_after) %yes

% Sanity Check
% Visualize r3d by itself (hist)
figure;histogram(temp.r3d,1000);xlabel('R_S_D');ylabel('N_m_e_a_s');
title('SD Separations (Matched2Data)');xlim([0 60])


%% Replace info.pairs with temp pairs and save pad file matched to data
info.pairs = temp;
save(['Pad_' info.optodes.CapName, '_matched2Data' '.mat'], 'info')


close all; clear all;

%% Appendix: Visualizing 2D layout in EEG style

% Load NIRx fullHead pad used in sparse pad light modeling tutorials
padname='FullHead_32x32';
load(['Pad_',padname,'.mat']);

% Turn eeg_style on and visualize
params_cap.eeg_style = 1;
PlotCap(info, params_cap)



