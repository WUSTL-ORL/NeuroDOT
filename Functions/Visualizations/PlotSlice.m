function PlotSlice(underlay, infoVol, params, overlay)

% PLOTSLICE Creates a plot of a 2D slice through 3D data. This function is
% essentially a sub-part of PlotSlices, so see that documentation for more
% information.
%
%
%   PLOTSLICE(underlay, infoVol) uses the volumetric data in "infoVol" to
%   display spatial coordinates of the slice in question.
%
%   PLOTSLICE(underlay, infoVol, params) allows the user to
%   specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       CH          1                       Turns crosshairs on (1) and off
%                                           (0).
%       Scale       (90% max)               Maximum value to which image is
%                                           scaled.
%       PD          0                       Declares that input image is
%                                           positive definite.
%       cbmode      0                       Specifies whether to use custom
%                                           colorbar axis labels.
%       cblabels    ([-90% max, 90% max])   Colorbar axis labels. When
%                                           cbmode==1, min defaults to 0 if
%                                           PD==1, both default to +/-
%                                           Scale if supplied. When
%                                           cbmode==0, then cblabels
%                                           dictates colorbar axis limits.
%       cbticks     (none)                  When cbmode==1, specifies
%                                           positions of tick marks on
%                                           colorbar axis.
%       slices      (center frames)         Select which slices are
%                                           displayed. If empty, activates
%                                           interactive navigation.
%       slices_type 'idx'                   Use MATLAB indexing ('idx') for
%                                           slices, or spatial coordinates
%                                           ('coord') as provided by
%                                           invoVol.
%       orientation 't'                     Select orientation of volume.
%                                           't' for transverse, 's' for
%                                           sagittal.
%       view                                which view of data to show; 
%                                           {'c','s','t'} for coronal,
%                                           sagittal, or transverse view.
%                                           Defaults to 't'.
%
%   Note: APPLYCMAP has further options for using "params" to specify
%   parameters for the fusion, scaling, and colormapping process.
%
%   PLOTSLICE(underlay, infoVol, params, overlay) overlays the image
%   provided by "overlay". When this is done, all color mapping is applied
%   to the overlay image, and the underlay is rendered as a grayscale image
%   times the RGB triplet in "params.BG".
%
% Dependencies: APPLYCMAP.
% 
% See Also: PLOTINTERPSURFMESH, PLOTSLICESMOV, PLOTSLICESTIMETRACE.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
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

%% Parameters and initialization.
LineColor = 'w';
BkgdColor = 'k';
[nVx, nVy, nVz] = size(underlay);
axlist = {'X', 'Y', 'Z'};

if ~exist('overlay', 'var')  ||  isempty(overlay)
    overlay = [];
else
    o_size = size(overlay);
    if (o_size(1) ~= nVx)  ||  (o_size(2) ~= nVy)  ||  (o_size(3) ~= nVz)
        error('*** Error: "underlay" size does not match "overlay" ***')
    end
end

if ~exist('infoVol','var')
    infoVol = [];
end

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [20, 200, 420, 420];
end
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end
if ~isfield(params, 'cboff')  
    params.cboff = 0;
end
if ~isfield(params, 'cbmode')  ||  isempty(params.cbmode)
    params.cbmode = 0;
end
% Set c_max. Ignore setting Scale, just ask if it's there.
if isfield(params, 'Scale')  &&  ~isempty(params.Scale)
    c_max = params.Scale;
elseif ~isempty(overlay)
    c_max = 0.9 * max(overlay(:));
else
    c_max = 0.9 * max(underlay(:));
end
% Set c_min and c_mid. Also ignore setting PD, which only matters for
% colormapping.
if isfield(params, 'PD')  &&  ~isempty(params.PD)  &&  params.PD
    c_mid = c_max / 2;
    c_min = 0;
else
    c_mid = 0;
    c_min = -c_max;
end
if params.cbmode
    if (~isfield(params, 'cbticks')  ||  isempty(params.cbticks))  &&...
            (~isfield(params, 'cblabels')  ||  isempty(params.cblabels))
        % If both are empty/not here, fill in defaults.
        params.cbticks = [0, 0.5, 1];
        params.cblabels = strtrim(cellstr(num2str([c_min, c_mid, c_max]', 3)));
    elseif ~isfield(params, 'cbticks')  ||  isempty(params.cbticks)
        % If only ticks missing...
        if numel(params.cblabels) > 2
            if isnumeric(params.cblabels)
                % If we have numbers, we can sort then scale them.
                params.cblabels = sort(params.cblabels);
                params.cbticks = (params.cblabels - params.cblabels(1))...
                    / (params.cblabels(end) - params.cblabels(1));
            else
                params.cbticks = 0:1/(numel(params.cblabels) - 1):1;
            end
        elseif numel(params.cblabels) == 2
            params.cbticks = [0, 1];
        else
            error('*** Need 2 or more colorbar ticks. ***')
        end
    elseif ~isfield(params, 'cblabels')  ||  isempty(params.cblabels)
        if numel(params.cbticks) > 2
            % If only labels missing, scale labels to tick spacing.
            scaled_ticks = (params.cbticks - params.cbticks(1)) /...
                (params.cbticks(end) - params.cbticks(1));
            params.cblabels = strtrim(cellstr(num2str(...
                [scaled_ticks * (c_max - c_min) + c_min]', 3)));
        elseif numel(params.cbticks) == 2
            params.cblabels = strtrim(cellstr(num2str([c_min, c_max]', 3)));
        else
            error('*** Need 2 or more colorbar labels. ***')
        end
    elseif numel(params.cbticks) == numel(params.cblabels)
        % As long as they match in size, continue on.
    else
        error('*** params.cbticks and params.cblabels do not match. ***')
    end
    % Any other errors are up to user.
else
    % Default cmin and cmax to work with "imagesc"; no cbticks needed.
    params.cblabels = [c_min, c_max];
end
if ~isfield(params, 'orientation')  ||  isempty(params.orientation)
    if isfield(infoVol, 'acq')
        params.orientation = infoVol.acq(1);
    else
        params.orientation = 't'; % This needs to come before params.slices!
    end
end
if ~isfield(params, 'slices')  ||  isempty(params.slices)
    S = [round(nVx / 2), round(nVy / 2), round(nVz / 2)];
else
    if isfield(params, 'slices_type')  &&  ~isempty(params.slices_type)
        switch params.slices_type
            case 'idx'
                % Do nothing.
            case 'coord'
                if exist('infoVol', 'var')  &&  isfield(infoVol, 'mmppix')
                    params.slices = change_space_coords(params.slices,...
                        infoVol, 'idx');
                else
                    error(['*** Error: No reference space provided ',...
                        'to calculate slice coordinates. ***'])
                end
        end
    end
    S = params.slices;
    % This corrects input for sagittal axis on 't' orientation.
    switch params.orientation
        case 't'
            S(1) = nVx - S(1) + 1;
    end
end
if strcmp(params.orientation, 't')
    underlay = flip(underlay, 1);
    if exist('overlay', 'var')
        overlay = flip(overlay, 1);
    end
end
if ~isfield(params, 'view'), params.view='t';end


%% Populate orientation stuff.
switch params.orientation
    case 's'
        orlist = {'Coronal', 'Transverse', 'Sagittal'};
        xdimlist = {nVz, nVz, nVx}; % derived from the volume, not a set of
        ydimlist = {nVy, nVx, nVy}; % canonical values.
        rot = {180, 180, 90};
        xlist = {'MR Dim 3' , 'MR Dim 3', 'MR Dim 1'};
        xflip = [1, 1, 0]; % DO NOT TOUCH! These flip the axes labels.
        ylist = {'MR Dim 2' , 'MR Dim 1', 'MR Dim 2'};
        yflip = [1, 1, 1];
    case 't'
        orlist = {'Sagittal', 'Coronal', 'Transverse'};
        xdimlist = {nVy, nVx, nVx};
        ydimlist = {nVz, nVz, nVy};
        rot = {90, 90, 90};
        xlist = {'MR Dim 2', 'MR Dim 1', 'MR Dim 1'};
        xflip = [0, 1, 1];
        ylist = {'MR Dim 3', 'MR Dim 3', 'MR Dim 2'};
        yflip = [1, 1, 1];
end

%% Apply color mapping.
if isempty(overlay)
    [FUSED, CMAP] = applycmap(underlay, [], params);
else
    [FUSED, CMAP] = applycmap(overlay, underlay, params);
end

%% Display the view

if exist('infoVol', 'var')  &&  isfield(infoVol, 'mmppix')
    coords = change_space_coords(S, infoVol, 'coord');
end

% Select view
switch params.view
    case 't'
        a=find(ismember(orlist,'Transverse'));        
    case 'c'
        a=find(ismember(orlist,'Coronal'));        
    case 's'
        a=find(ismember(orlist,'Sagittal'));
end

% Acquire overlay image.
switch a
    case 1
        SLICE = squeeze(FUSED(S(a), :, :, :));
    case 2
        SLICE = squeeze(FUSED(:, S(a), :, :));
    case 3
        SLICE = squeeze(FUSED(:, :, S(a), :));
end

% Display image.
if params.cbmode
    image(imrotate(SLICE, rot{a}))
else
    imagesc(imrotate(SLICE, rot{a}), params.cblabels)
end

% Titles and labels.
axis image
h = gca;
set(h, 'XColor', LineColor, 'YColor', LineColor, 'Box', 'on');

if strcmp(params.orientation, 's')  &&  ~strcmp(orlist{a}, 'Sagittal')
    cel = {[orlist{a}, ' View']; 'Left is Left'; ['Frame ', num2str(S(a))]};
elseif strcmp(params.orientation, 't')  &&  strcmp(orlist{a}, 'Sagittal')
    cel = {[orlist{a}, ' View']; ['Frame ', num2str(nVx - S(a) + 1)]};
else
    cel = {[orlist{a}, ' View']; ['Frame ', num2str(S(a))]};
end

if exist('coords', 'var')
    if strcmp(params.orientation, 's')  &&  ~strcmp(orlist{a}, 'Sagittal')
        cel{end + 1} = [axlist{a} ' = ' num2str(coords(a))];
    elseif strcmp(params.orientation, 't')  &&  strcmp(orlist{a}, 'Sagittal')
        cel{end + 1} = [axlist{a} ' = ' num2str(coords(a).*-1)];
    else
        cel{end + 1} = [axlist{a} ' = ' num2str(coords(a))];
    end
end

title(cel, 'Color', LineColor, 'FontSize', 12)
xlabel(xlist{a}, 'Color', LineColor, 'FontSize', 10)
ylabel(ylist{a}, 'Color', LineColor, 'FontSize', 10)

% Flip axes tick labels as necessary.
if  xflip(a)
    h.XTickLabel = flip(h.XTickLabel); % flips the labels but not the data.
    incr = mean(diff(h.XTick));
    dif = incr + max(h.XTick(:)) - xdimlist{a}; % fixes the positions of the labels.
    h.XTick = h.XTick - dif  + 1;
    h.XTickLabel(h.XTick == 1) = [];
    h.XTick(h.XTick == 1) = [];
end
if  yflip(a)
    h.YTickLabel = flip(h.YTickLabel);
    incr = mean(diff(h.YTick));
    dif = incr + max(h.YTick(:)) - ydimlist{a};
    h.YTick = h.YTick - dif + 1;
    h.YTickLabel(h.YTick == 1) = [];
    h.YTick(h.YTick == 1) = [];
end

% Add a colorbar.
colormap(gca,CMAP)
h2 = colorbar('Color', LineColor);
if params.cbmode
    set(h2, 'Ticks', params.cbticks, 'TickLabels', params.cblabels);
end
if params.cboff
   set(h2,'Visible','off') 
end

%
