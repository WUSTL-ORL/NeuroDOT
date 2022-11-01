%% Program Description
%
% This script provides some examples to illustrate how to use the
% PlotCapMeasSpaceData_Pairs function and some of its features.
% 
% Copyright (c) 2017 Washington University 
% Created By: Zachary E. Markow, Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
%


%% Load a cap pad file and create some example data.

%padFile = load('Pad_Adult_96x92.mat');
padFile = load('Pad_Adult_96x92_on_Example_Mesh.mat');

dimString = '2D';  % '2D' or '3D' depending on pad file and what user wants to display.

nMeasTotal = length(padFile.info.pairs.Src);

randnData = randn(nMeasTotal,1);

% Pretend that all measurements involving detectors 76-80 are the only bad
% measurements.
padFile.info.MEAS = [];
padFile.info.MEAS.GI = ~ismember(padFile.info.pairs.Det, 76:80);

% Pretend that all measurements involving sources 37-41 have value = 3.
randnData(ismember(padFile.info.pairs.Src,37:41)) = 3;

% Pretend that all measurements involving sources 57-61 have value = -2.
randnData(ismember(padFile.info.pairs.Src,57:61)) = -2;

randnDataAbs = abs(randnData);


%% Set color bar thresholds = 0, maximum value = 4, and only display measurements with 2D source-detector distance < 50 mm.

% By default, only good measurements will be displayed.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Repeat, but set thresholds and custom tick labels on colormap.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.Th.P = 1.0;
params.Th.N = -0.5;
params.Scale = 4;
params.cbticks_in_data_units = [-4 -3 -2 -0.5 0 1.0 1.5 2 3.5 4];
params.rlimits = [0,50];
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Repeat but with zero thresholds and show only bad measurements within the selected source-detector distance range.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];
params.mode = 'bad';
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Add custom title and do not mark bad optodes with circles.

params = [];  % Reset params.
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.title = 'Example Data from randn Function';
params.markBadOptodes = 0;
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];


%% Use custom colorbar tick marks and a custom colormap.  Also display good and bad measurments.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.cbticks_in_data_units = [-3.5 -2 0 1 3 4];
params.Cmap = 'cool';
params.mode = 'both';
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Use custom colorbar tick marks and text labels.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.cbticks_in_data_units = [-3.5 -0.5 3.0];
params.cblabels = {'very low', 'near mid', 'high'};
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Display with colorbar to the left of the cap.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.CBar_location = 'WestOutside';  % Can be 'SouthOutside' to place below cap.
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,50];
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Turn off colorbar and custom title, restrict to measurements with shorter source-detector separation, and specify line thickness.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.CBar_on = 0;
params.title = '';
params.Th.P = 0;
params.Scale = 4;
params.rlimits = [0,30];
params.line_thickness = 2.0;
PlotCapMeasSpaceData_Pairs(randnData, padFile.info, params);


%% Show good and bad measurements with green and red lines.

% Note that to display a line for both types of measurements, params.mode
% must be set to 'both', but this causes the title's reported percentage of
% measurements that are shown in the specified source-detector range to be
% 100%.  This differs from the behavior of PlotCapGoodMeas.

params = [];
params.dimension = dimString;  % Pull dimString from above so that we do not have to make new copies of this script to work with 3D case.
params.mode = 'both';
params.Cmap = [1.0 0 0; 0 1.0 0];
params.DR = size(params.Cmap,1);
params.PD = 1;
params.Scale = 1;
params.Th.P = 0;
params.line_thickness = 1.0;
params.rlimits = [0,50];
params.cbticks_in_data_units = sort([0:0.25:1.0, 0.35, 0.85]);
goodBad = repmat(0.01, nMeasTotal, 1);  % Assume bad initially.
goodBad(padFile.info.MEAS.GI) = 1.0;  % Replace good with values of 1.0.
PlotCapMeasSpaceData_Pairs(goodBad, padFile.info, params);


%% Check that wavelength averaging is working.

nonrandData = repmat(0.5,nMeasTotal,1);
nonrandData(padFile.info.pairs.WL == 2) = 1.5;
infoTemp = padFile.info;
infoTemp.MEAS.GI = infoTemp.MEAS.GI & ~(infoTemp.pairs.WL == 2 & infoTemp.pairs.Src == 10);
params = [];
%params.mode = 'good';
params.Cmap = 'jet';
params.PD = 1;
params.Scale = 2;
params.Th.P = 0;
params.rlimits = [0,40];
params.cbticks_in_data_units = 0:0.5:2;
PlotCapMeasSpaceData_Pairs(nonrandData, infoTemp, params);

% Try 3 subplots, one with each wavelength, and one with averaging over
% wavelengths.
figure('Color','k');
h1 = subplot(3,1,1);
params.WLs = 1;
params.fig_handle = h1;
PlotCapMeasSpaceData_Pairs(nonrandData, infoTemp, params);
h2 = subplot(3,1,2);
params.WLs = 2;
params.fig_handle = h2;
PlotCapMeasSpaceData_Pairs(nonrandData, infoTemp, params);
h3 = subplot(3,1,3);
params.WLs = [1,2];
params.fig_handle = h3;
PlotCapMeasSpaceData_Pairs(nonrandData, infoTemp, params);
