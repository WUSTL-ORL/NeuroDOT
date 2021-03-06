%% Script for generating a Pad file

% The Pad file contains information about your source-detector array: where
% the optical elements are in space, the topology of the grid, information
% about the measurements, and information about the tissue model, if
% applicable. All of this information is contained within a structure of a
% variable called info.
% This script will walk through how to generate a Pad file for
% efficient data processing in NeuroDOT. 

%% info.optodes

% If you already have an array design, populate the information as below.
% In this example, grid is a structure containing 2D and 3D versions of the
% source/detector locations.
info.optodes.CapName='Your_Cap_Name';
info.optodes.dpos2=grid.dpos2; % 2D detector locations
info.optodes.dpos3=grid.dpos3; % 3D detector locations
info.optodes.spos2=grid.spos2; % 2D source locations
info.optodes.spos3=grid.spos3; % 3D source locations

% Define the following as given; this will be updated in a future release
info.optodes.plot3orientation.i='R2L';  
info.optodes.plot3orientation.j='P2A'; 
info.optodes.plot3orientation.k='D2V'; 

% If you need to generate the topology structure for your grid:
dr=5;       % Distance separating 'nearest neighbor' groups
Rad=Grid2Radius_180824(grid,dr);


%% If instead you are designing a grid from scratch in style of our visual cap:
Ndy=4;          % Num detectors in y-direction
Ndx=7;          % Num detectors in x-direction
NN1sep=13;      % Closest source-detector Separation
Srad=80;        % Radius of spherical head model
[grid,Rad]=Make_Sphere_Cap(Ndy,Ndx,NN1sep,Srad);

% grid contains information on the positions of the sources and detectors.
% Rad contains information on the topology of the array.

PlotSD(grid.spos2,grid.dpos2,'norm'); % Visualize 2D positions
PlotSD(grid.spos3,grid.dpos3,'norm'); % Visualize 3D positions



%% info.pairs

% This part of the structure contains the full measurement list for your
% array in a table with the following columns (each row is a measurement):
% *.Src     source number
% *.Det     detector number
% *.NN      nearest neighbor type
% *.WL      wavelength number index
% *.lambda  actual wavelength for measurement
% *.Mod     modulation type; 'CW' is for continuous wave;
% *.r2d     2-dimentional source-detector separation
% *.r3d     3-dimentional source-detector separation

% Once you have the grid and radius set up (as above), this structure can
% be generated with the following:
grid.name='Your_Grid_Name';
params.lambda=[750,850];    % Wavelengths for your system
params.Mod='CW';            % Modulation strategy
[info]=Generate_Info_from_Grid_Rad(grid,Rad,params);



%% info.tissue

% This structure contains information about the light model you will
% generate. 
info.tissue.affine=eye(4); % if affine transform to atlas space is Identity
info.tissue.affine_target='MNI'; % Atlas target from your model

% These next pieces will be created when you make the light model
% info.tissue.info=infoT1; % Spatial metadata for a head model, on which the grid will be relaxed
% info.tissue.dim % spatial metadata for sub-space containing Jacobian



%% Then just navigate to Support_Files/Pad and save your pad
tag='_180901';
save(['Pad_',info.optodes.CapName,tag],'info')







