function PlotSlices(underlay, infoVol, params, overlay)

% PLOTSLICES Creates an interactive three-slice plot.
%
%   PLOTSLICES(underlay) takes a 3D voxel space image "underlay" and
%   generates views along the three canonical axes.
%
%   In interactive mode, left-click on any point to move to those slices.
%   To reset to the middle of the volume, right-click anywhere. To cancel
%   interactive mode, press "Q", "Esc", or the middle mouse button.
%
%   PLOTSLICES(underlay, infoVol) uses the volumetric data in "infoVol" to
%   display spatial coordinates of the slices in question.
%
%   PLOTSLICES(underlay, infoVol, params) allows the user to
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
%       cboff       0                       If set to 1, no colorbar is
%                                           displayed
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
%       bgW        [0,0,0]                  For images with an overlay, change
%                                           the background color.
%                                           Options are all RGB values [0,1].
%                                           For paper-ready figures with a
%                                           white background [1,1,1].
%
%   Note: APPLYCMAP has further options for using "params" to specify
%   parameters for the fusion, scaling, and colormapping process.
%
%   PLOTSLICES(underlay, infoVol, params, overlay) overlays the image
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
button = 0;
axlist = {'X', 'Y', 'Z'};

% Set simple default for single input case
if nargin<3
    params.PD=1;
    params.Cmap.P=gray(1000);
    params.Th.P=0;
    params.Th.N=0;
    params.Scale=max(underlay(:));
    params.BG=[0,0,0];
end

if ~exist('params', 'var')
    params = [];
end

% if ~isfield(params, 'PD'),      params.PD=1;end
% if ~isfield(params, 'Cmap'),    params.Cmap.P=gray(1000);end
% if ~isfield(params, 'Th'),      params.Th.P=0;params.Th.N=0;end
% if ~isfield(params, 'Scale'),   params.Scale=max(underlay(:));end
% if ~isfield(params, 'BG'),      params.BG=[0,0,0];end

if ~exist('infoVol','var')
    infoVol = [];
end

if ~exist('overlay', 'var')  ||  isempty(overlay)
    overlay = [];
else
    o_size = size(overlay);
    if (o_size(1) ~= nVx)  ||  (o_size(2) ~= nVy)  ||  (o_size(3) ~= nVz)
        error('*** Error: "underlay" size does not match "overlay" ***')
    end
    if isfield(params,'BG'),params=rmfield(params,'BG');end
end


if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [20, 200, 1240, 420];
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
if ~isfield(params, 'CH')  ||  isempty(params.CH)
    params.CH = 1;
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
            params.cblabels = strtrim(cellstr(num2str([scaled_ticks * (c_max - c_min) + c_min]', 3)));
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
% Adjust aspect ratio relative to dim 1
if isfield(infoVol,'mmx')
    arX=infoVol.mmx./infoVol.mmx;
    arY=infoVol.mmx./infoVol.mmy;
    arZ=infoVol.mmx./infoVol.mmz;
else
    arX=1;
    arY=1;
    arZ=1;
end
if ~isfield(params, 'slices')  ||  isempty(params.slices)
    S = [round(nVx / 2), round(nVy / 2), round(nVz / 2)];
    lookon = 1;
else
    if isfield(params, 'slices_type')  &&  ~isempty(params.slices_type)
        switch params.slices_type
            case 'idx'
                % Do nothing.
            case 'coord'
                if exist('infoVol', 'var')  &&  isfield(infoVol, 'mmppix')
                    params.slices = change_space_coords(params.slices, infoVol, 'idx');
                else
                    error('*** Error: No reference space provided to calculate slice coordinates. ***')
                end
        end
    end
    S = params.slices;
    % This corrects input for sagittal axis on 't' orientation.
    switch params.orientation
        case 't'
            S(1) = nVx - S(1) + 1;
    end
    lookon = 0;
end
if strcmp(params.orientation, 't')
    underlay = flip(underlay, 1);
    if exist('overlay', 'var')
        overlay = flip(overlay, 1);
    end
end

if ~isfield(params, 'bgW')
    params.bgW = [0,0,0];
end

%% Populate orientation stuff.
switch params.orientation
    case 's'
        orlist = {'Coronal', 'Transverse', 'Sagittal'};
        dimlist = {nVx, nVy, nVz}; % Note, this corresponds to the dimensions
        xdimlist = {nVz, nVz, nVx}; % derived from the volume, not a set of
        ydimlist = {nVy, nVx, nVy}; % canonical values.
        rot = {180, 180, 90};
        CHs = {'plot(ones(1, nVy) .* (nVz - S(3) + 1), 1:nVy, ''-m'', ''LineWidth'', 1)',...
            'plot(1:nVz, ones(1, nVz) .* (nVy - S(2) + 1), ''-g'', ''LineWidth'', 1)';...
            'plot(ones(1, nVx) .* (nVz - S(3) + 1), 1:nVx, ''-m'', ''LineWidth'', 1)',...
            'plot(1:nVz, ones(1, nVz) .* (nVx - S(1) + 1), ''-c'', ''LineWidth'', 1)';...
            'plot(ones(1, nVy) .* S(1), 1:nVy, ''-c'', ''LineWidth'', 1)',...
            'plot(1:nVx, ones(1, nVx) .* (nVy - S(2) + 1), ''-g'', ''LineWidth'', 1)'};
        xlist = {'MR Dim 3' , 'MR Dim 3', 'MR Dim 1'};
        xflip = [1, 1, 0]; % DO NOT TOUCH! These flip the axes labels.
        ylist = {'MR Dim 2' , 'MR Dim 1', 'MR Dim 2'};
        yflip = [1, 1, 1];
        arXYZ=[arZ,arY,arX;arZ,arX,arY;arX,arY,arZ];
        
    case 't'
        orlist = {'Sagittal', 'Coronal', 'Transverse'};
        dimlist = {nVx, nVy, nVz};
        xdimlist = {nVy, nVx, nVx};
        ydimlist = {nVz, nVz, nVy};
        rot = {90, 90, 90};
        CHs = {'plot(ones(1, nVz) .* S(2), 1:nVz, ''-c'', ''LineWidth'', 1)',...
            'plot(1:nVy, ones(1, nVy) .* (nVz - S(3) + 1), ''-g'', ''LineWidth'', 1)';...
            'plot(ones(1, nVz) .* S(1), 1:nVz, ''-m'', ''LineWidth'', 1)',...
            'plot(1:nVx, ones(1, nVx) .* (nVz - S(3) + 1), ''-g'', ''LineWidth'', 1)';...
            'plot(ones(1, nVy) .* S(1), 1:nVy, ''-m'', ''LineWidth'', 1)',...
            'plot(1:nVx, ones(1, nVx) .* (nVy - S(2) + 1), ''-c'', ''LineWidth'', 1)'};
        xlist = {'MR Dim 2', 'MR Dim 1', 'MR Dim 1'};
        xflip = [0, 1, 1];
        ylist = {'MR Dim 3', 'MR Dim 3', 'MR Dim 2'};
        yflip = [1, 1, 1];
        arXYZ=[arY,arZ,arX;arX,arZ,arY;arX,arY,arZ];
end

%% Apply color mapping.
if isempty(overlay)
    [FUSED, CMAP] = applycmap(underlay, [], params);
    if isfield(params, 'bgW')
        FUSED = reshape(FUSED, [], 3);
        idx = sum(~FUSED,2)==3;
        FUSED(idx,:) = repmat([params.bgW(1),params.bgW(2),params.bgW(3)], sum(idx), 1);
        FUSED = reshape(FUSED, [infoVol.nVx, infoVol.nVy, infoVol.nVz,3]);
    end
    if isempty(FUSED),return;end
    
else
    [FUSED, CMAP] = applycmap(overlay, underlay, params);
    if isfield(params, 'bgW')
        FUSED = reshape(FUSED, [], 3);
        idx = sum(~FUSED,2)==3;
        FUSED(idx,:) = repmat([params.bgW(1),params.bgW(2),params.bgW(3)], sum(idx), 1);
        FUSED = reshape(FUSED, [infoVol.nVx, infoVol.nVy, infoVol.nVz,3]);
    end
    if isempty(FUSED),return;end
    
end



%% Display the views on three subplots in a while loop for point-and-click navigation.
N = 0;
disp(['PlotSlices is active - Please press q or the middle mouse button to quit'])
while ~any(button == [2, 27, 81, 113]) % 2 = middle mouse button, 27 = Esc, 81 = Q, 113 = q.
    try
        N = N + 1;
        if exist('infoVol', 'var')  &&  isfield(infoVol, 'mmppix')
            coords = change_space_coords(S, infoVol, 'coord');
        end
        for ax = {1, 2, 3}
            a = ax{1};
            
            % Create subplot.
            subplot(1, 3, a)
            
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
            if (N == 1)  &&  xflip(a)
                h.XTickLabel = flip(h.XTickLabel); % This is mainly window dressing, flips the labels but not the data.
                incr = mean(diff(h.XTick));
                dif = incr + max(h.XTick(:)) - xdimlist{a}; % As such, this formula fixes the positions of the labels.
                h.XTick = h.XTick - dif  + 1;
                h.XTickLabel(h.XTick == 1) = [];
                h.XTick(h.XTick == 1) = [];
            end
            if (N == 1)  &&  yflip(a)
                h.YTickLabel = flip(h.YTickLabel);
                incr = mean(diff(h.YTick));
                dif = incr + max(h.YTick(:)) - ydimlist{a};
                h.YTick = h.YTick - dif + 1;
                h.YTickLabel(h.YTick == 1) = [];
                h.YTick(h.YTick == 1) = [];
            end
            
            % Fix Aspect Ratio
            daspect(arXYZ(a,:))
            
            % This saves the axes as an object to be used later in
            % navigation.
            eval([lower(orlist{a}(1:3)), 'ax = gca;'])
            
            % Add crosshairs.
            if params.CH
                hold on
                eval(CHs{a, 1})
                eval(CHs{a, 2})
            end
        end
        
       % Add a colorbar.
        colormap(CMAP)
        h2 = colorbar(eval([lower(orlist{end}(1:3)), 'ax']), 'Color', LineColor);
        
        % Py version has section here where we double check that ticks are kosher
        % Not sure if this is 100% necessary for matlab
        ticks = h2.Ticks; %get whatever ticks matlab is using
        ticks_new  = [min(ticks), ((max(ticks)-min(ticks))/2) + min(ticks), max(ticks)];
        if ticks_new(2) > ticks_new(3) %if mid > max, reverse sign of mid
            ticks_new(2) = -1 * ticks_new(2);
        end
        limits = h2.Limits;
        if ticks_new(1) < limits(1) %if lowest tick less than lower limit of cb, set lowest tick to lower limit of cb
            ticks_new(1) = limits(1);
        end
        if ticks_new(3) > limits(2) %if highest tick greater than upper limit of cb, set highest tick to upper limit of cb
            ticks_new(3) = limits(2);
        end
        set(h2, 'Ticks', ticks_new);
        % Set labels to min, mid, and max of colorscale
        set(h2, 'TickLabels', {round(c_min,3); round(c_mid,4); round(c_max,3)})
        
        % Back to OG colorbar code
        if params.cbmode
            set(h2, 'Ticks', params.cbticks, 'TickLabels', params.cblabels);
        end
        if params.cboff
           set(h2,'Visible','off') 
        end
        
        %% Add point-and-click navigation.
        if lookon == 0
            break
        end
        oldS = S;
        [gx, gy, button] = ginput(1);
        gx = round(gx);
        gy = round(gy);
        if gx < 1, gx = 1;end
        if gy < 1, gy = 1;end
        if button == 1
            switch params.orientation
                case 't'
                    switch gca
                        case sagax
                            if gx > nVy, gx = nVy;end
                            if gy > nVz, gy = nVz;end
                            S(2) = gx;
                            S(3) = nVz - gy + 1;
                        case corax
                            if gx > nVx, gx = nVx;end
                            if gy > nVz, gy = nVz;end
                            S(1) = gx;
                            S(3) = nVz - gy + 1;
                        case traax
                            if gx > nVx, gx = nVx;end
                            if gy > nVy, gy = nVy;end
                            S(1) = gx;
                            S(2) = nVy - gy + 1;
                    end
                case 's'
                    switch gca
                        case sagax
                            if gx > nVx, gx = nVx;end
                            if gy > nVy, gy = nVy;end
                            S(1) = gx;
                            S(2) = nVy - gy + 1;
                        case corax
                            if gx > nVz, gx = nVz;end
                            if gy > nVy, gy = nVy;end
                            S(2) = nVy - gy + 1;
                            S(3) = nVz - gx + 1;
                        case traax
                            if gy > nVx, gy = nVx;end
                            if gx > nVz, gx = nVz;end
                            S(1) = nVx - gy + 1;
                            S(3) = nVz - gx + 1;
                    end
            end
        elseif button == 3
            S(1) = round(nVx / 2);
            S(2) = round(nVy / 2);
            S(3) = round(nVz / 2);
        end
    catch err
%         rethrow(err)
        S = oldS;
    end
end



%
