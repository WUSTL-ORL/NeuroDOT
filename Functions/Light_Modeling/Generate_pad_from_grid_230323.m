function info = Generate_pad_from_grid_230323(grid, params, info)

% This function generates info structures optodes and pairs from the 
% data in grid and radius files. 
% First a radius structure is created and populated with topological 
% information for a cap.

% The input: "grid" must have fields that list the spatial locations of 
% sources and detectors in 3D: spos3, dpos3.

% The input: "params" can be used to pass in mod type (default is 'CW' 
% but can be the modulation frequency if fd) and wavelength(s) of the 
% data in field 'lambda'.
%
% The params input can be used to populate the following information:
%   params.dr           Minimum separation for sources and detectors 
%                       to be grouped into different neighbors.
%   params.lambda       Wavelengths of the light. Any number of comma-
%                       separated values is allowed. Default: [750,850].
%   params.Mod          Modulation type or frequency. Can be 'CW or 'FD' or
%                       can be the actual modulation frequency (e.g., 0 or
%                       200) in MHz.
%   params.pos2         Indicator. Determines NN classification. Defaults
%                       to 0, where 3D coordinates are used. If set to 1,
%                       2D coordinates will be used for NN classification.
%   params.CapName      Name for your pad file.

% The input: "info" is optional, as this function will create a info
% structure. If an info structure is passed in as an input, the
% sub-structures info.optodes and info.pairs will be overwritten. All
% information outside of info.optodes and info.pairs will be preserved.


%% Parameters and Initialization

% Params and output structure
if ~exist('info', 'var'), info = struct;end
if ~exist('params', 'var'), params = struct;end
if ~isfield(params, 'dr'),dr=10; % defaults to 10mm
else
    dr = params.dr;
end
if ~isfield(params, 'lambda'), params.lambda = [750,850];end
if ~isfield(params, 'Mod'), params.Mod = 'CW';end
if ~isfield(params, 'pos2'), params.pos2 = 0;end

% Optode positions
% Make generalizable (can be either 2D or 3D) optode pos
% Defaults to 3D pos if 2D pos not part of input grid structure
if ~isfield(grid, 'spos')
    if isfield(grid, 'spos2')
        grid.spos = grid.spos2;
    else
        grid.spos = grid.spos3;
    end
end
if ~isfield(grid, 'dpos')
    if isfield(grid, 'dpos2')
        grid.dpos = grid.dpos2;
    else
        grid.dpos = grid.dpos3;
    end
end
% If 3D not supplied as input, set to 2D pos where col3 is all zeros
if ~isfield(grid, 'spos3')
    grid.spos3 = cat(2,grid.spos,zeros(size(grid.spos,1),1));
end
if ~isfield(grid, 'dpos3')
    grid.dpos3 = cat(2,grid.dpos,zeros(size(grid.dpos,1),1));
end

% Calculate number of sources, detectors, measurements, and wavelengths
Ns = size(grid.spos3,1);
Nd = size(grid.dpos3,1);
Nm = Ns*Nd;
Nwl = length(params.lambda);

% Initialize SD separations and measurement list
% Note: these are all for a single wavelength, and will get replicated 
%   in a below section of the function for all other WL
r2 = zeros(Nm,1); %2D SD separations
r3 = zeros(Nm,1); %3D SD separations
measList = zeros(Nm,2); %basic measurement list that only has [Src, Det]


%% Populate info.optodes structure

% Detectors
if isfield(params,'CapName'),info.optodes.CapName = params.CapName;end
if isfield(grid,'dpos3'),info.optodes.dpos3 = grid.dpos3;end
if isfield(grid,'dpos2')
    info.optodes.dpos2 = grid.dpos2;
elseif isfield(grid,'dpos')
    info.optodes.dpos2 = grid.dpos;
else
    info.optodes.dpos2 = grid.dpos3;
end

% Sources
if isfield(grid,'spos3'),info.optodes.spos3 = grid.spos3;end
if isfield(grid,'spos2')
    info.optodes.spos2 = grid.spos2;
elseif isfield(grid,'spos')
    info.optodes.spos2 = grid.spos;
else
    info.optodes.spos2 = grid.spos3;
end


%% Make Measlist, r3d, and r2d

m = 0;
for d = 1:Nd
    for s = 1:Ns
        m = m+1;
        measList(m,1) = s;
        measList(m,2) = d;
        r2(m) = norm(info.optodes.spos2(s,:)-info.optodes.dpos2(d,:));
        r3(m) = norm(info.optodes.spos3(s,:)-info.optodes.dpos3(d,:));
        
    end
end


%% Populate info.pairs structure

info.pairs = struct;
info.pairs.Src = repmat(measList(:,1),Nwl,1);
info.pairs.Det = repmat(measList(:,2),Nwl,1);
% info.pairs.NN will be created and populated below
info.pairs.WL = nan(size(info.pairs.Src));
info.pairs.lambda = nan(size(info.pairs.Src));
info.pairs.Mod = repmat(params.Mod,Nm*Nwl,1);
info.pairs.r2d = repmat(r2,Nwl,1);
info.pairs.r3d = repmat(r3,Nwl,1);

% Make sure WL and Lambda reflect both wavelengths
for j = 1:Nwl
    info.pairs.WL(Nm*(j-1)+1:Nm*j) = ones(Nm,1).*j;
    info.pairs.lambda(Nm*(j-1)+1:Nm*j) = params.lambda(j);
end


%% Populate NN's 

info = calc_NN(info,params.dr);


end