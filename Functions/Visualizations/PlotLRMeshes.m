function [hPatchesL,hPatchesR,params2,hF] = PlotLRMeshes(meshL, meshR, params)
% PLOTLRMESHES Renders a pair of hemispheric meshes.
%
%   PLOTLRMESHES(meshL, meshR) renders the data in a pair of left and right
%   hemispheric meshes "meshL.data" and "meshR.data", respectively, and
%   applies full color mapping to them. If no data is present for either
%   mesh, a default gray mesh will be plotted.
%
%   PLOTLRMESHES(meshL, meshR, params) allows the user to
%   specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 960, 420]     Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       Scale       (90% max)               Maximum value to which image is
%                                           scaled.
%       PD          0                       Declares that input image is
%                                           positive definite.
%       cblabels    ([-90% max, 90% max])   Colorbar axis labels. Min
%                                           defaults to 0 if PD==1, both
%                                           default to +/- Scale if
%                                           supplied.
%       cbticks     (none)                  Specifies positions of tick
%                                           marks on colorbar axis.
%       alpha       1                       Transparency of mesh.
%       view        'lat'                   Sets the view perspective.
%
% Dependencies: ADJUST_BRAIN_POS, ROTATE_CAP, ROTATION_MATRIX, APPLYCMAP.
%
% See Also: PLOTINTERPSURFMESH, VOL2SURF_MESH.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Other Contributor(s): Zachary E. Markow
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

%% Parameters and Initialization.
LineColor = 'w';
BkgdColor = 'k';
FaceColor = 'interp';

if ~isfield(meshL, 'data')  ||  isempty(meshL.data)
    meshL.data = zeros(size(meshL.nodes, 1), 1);
end
if ~isfield(meshR, 'data')  ||  isempty(meshR.data)
    meshR.data = zeros(size(meshR.nodes, 1), 1);
end

if (size(meshL.data, 2) > 1)  ||  (size(meshR.data, 2) > 1)
    error(['*** One or both mesh imputs has more than one time points in',...
        ' "mesh.data". PlotLRMeshes can only plot one time point. ***'])
end

if ~exist('params', 'var')  ||  isempty(params)
    params = [];
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [20, 200, 960, 420];
end
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            parentFig = get(params.fig_handle,'Parent');
            set(groot, 'CurrentFigure', parentFig);
            set(parentFig, 'CurrentAxes', params.fig_handle);
    end
end
% Slotted this in before 'Scale' to create the "gray mode" when no data is
% input.
if ~any([meshL.data(:); meshR.data(:)])
    params.Scale = 1;
end
if ~isfield(params, 'Scale')  ||  isempty(params.Scale)
    params.Scale = 0.9 * max([meshL.data(:); meshR.data(:)]);
end
if ~isfield(params, 'Th')  ||  isempty(params.Th)
    params.Th.P = 0.25 * params.Scale;
    params.Th.N = -params.Th.P;
end
if ~isfield(params, 'alpha')  ||  isempty(params.alpha)
    params.alpha = 1;
end
if ~isfield(params, 'lighting')  ||  isempty(params.lighting)
    params.lighting = 'gouraud';
end
if ~isfield(params, 'ctx')  ||  isempty(params.ctx)
    params.ctx = 'std';
end
if ~isfield(params, 'view')  ||  isempty(params.view)
    params.view = 'lat';
    if strcmp(params.ctx, 'flat')
        params.view = 'lat';
    end
end

if ~isfield(params, 'FaceColor')
    params.FaceColor = [0.5,0.5,0.5];
end
% Set c_max. Ignore setting Scale, just ask if it's there.
if isfield(params, 'Scale')  &&  ~isempty(params.Scale)
    c_max = params.Scale;
else
    c_max = 0.9 * max([meshL.data(:); meshR.data(:)]);
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

if ~isfield(params,'CBar_on'),params.CBar_on=0;end 
if params.CBar_on
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
            params.cblabels = ...
                strtrim(cellstr(num2str(...
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
    end % Any other errors are up to user.
end

%% Image Hemispheres. Reposition meshes for standard transverse orientation.
[Lnodes, Rnodes] = adjust_brain_pos(meshL, meshR, params);


%% Set Lighting and Persepctive.
if isfield(params, 'render_tool')
    delete(gca);
end
switch params.view
    case {'lat', 'med'}
        view([-90, 0]);
        light('Position', [-100, 200, 0], 'Style', 'local');
        % These two lines create minimal lighting good luck.
        light('Position', [-50, -500, 100], 'Style', 'infinite');
        light('Position', [-50, 0, 0], 'Style', 'infinite');
    case {'post', 'dorsal'} % [x:L -> R, y:P -> A, z:V -> D]
        if strcmp(params.view, 'post')
            view([0 0]);
        end
        if strcmp(params.view, 'dorsal')
            view([0 90]);
            light('Position', [100, 300, 100], 'Style', 'infinite');
        end
        light('Position', [-500, -20, 0], 'Style', 'local');
        light('Position', [500, -20, 0], 'Style', 'local');
        light('Position', [0, -200, 50], 'Style', 'local');
end
if isfield(params, 'render_tool')
    set(params.fig_handle, 'Units', 'pixels');
end


%% Adjust visualization parameters
switch params.ctx
    case 'flat'
        params.FaceColor = 'interp';
        FaceColor = params.FaceColor;
        params.EdgeColor = 'interp';
        EdgeColor = params.EdgeColor;
        if isfield(params, 'reg')
            if params.reg
                [FV_CDataL,CMAP,params2] = applycmap(meshL.data, meshL.region, params);
                [FV_CDataR,CMAP,params2] = applycmap(meshR.data, meshR.region, params);
            else
                params.FaceColor = 'interp'; params.EdgeColor = 'None';
                [FV_CDataL, CMAP, params2] = applycmap(meshL.data, [], params);
                [FV_CDataR, CMAP, params2] = applycmap(meshR.data, [], params);
            end
        else
            params.FaceColor = 'interp'; params.EdgeColor = 'None';
                [FV_CDataL, CMAP, params2] = applycmap(meshL.data, [], params);
                [FV_CDataR, CMAP, params2] = applycmap(meshR.data, [], params);
        end
            hPatchesL = patch('Faces', meshL.elements(:, 1:3), 'Vertices', Lnodes,...
            'EdgeColor', EdgeColor, 'FaceColor', FaceColor, 'FaceVertexCData', FV_CDataL,...
            'FaceLighting', params.lighting, 'FaceAlpha', params.alpha,...
            'AmbientStrength', 0.25, 'DiffuseStrength', .75, 'SpecularStrength', .1);

        % Image Right Side.
        hPatchesR = patch('Faces', meshR.elements(:, 1:3), 'Vertices', Rnodes,...
            'EdgeColor', EdgeColor, 'FaceColor', FaceColor, 'FaceVertexCData', FV_CDataR,...
            'FaceLighting', params.lighting, 'FaceAlpha', params.alpha,...
            'AmbientStrength', 0.25, 'DiffuseStrength', .75, 'SpecularStrength', .1);

        
    case {'std', 'inf', 'vinf'}
        %% Image Left Side.
        [dataL, CMAP, params2] = applycmap(meshL.data, [], params);
        hPatchesL = patch('Faces', meshL.elements(:, 1:3), 'Vertices', Lnodes,...
            'EdgeColor', 'none', 'FaceColor', FaceColor, 'FaceVertexCData', dataL,...
            'FaceLighting', params.lighting, 'FaceAlpha', params.alpha,...
            'AmbientStrength', 0.25, 'DiffuseStrength', .75, 'SpecularStrength', .1);

        % Image Right Side.
        [dataR, CMAP, params2] = applycmap(meshR.data, [], params);
        hPatchesR = patch('Faces', meshR.elements(:, 1:3), 'Vertices', Rnodes,...
            'EdgeColor', 'none', 'FaceColor', FaceColor, 'FaceVertexCData', dataR,...
            'FaceLighting', params.lighting, 'FaceAlpha', params.alpha,...
            'AmbientStrength', 0.25, 'DiffuseStrength', .75, 'SpecularStrength', .1);

end



%% Visual Formatting.
set(gca, 'Color', BkgdColor, 'XColor', LineColor, 'YColor', LineColor);

axis image
axis off
% 
% title('Volumetric Surface Mapping', 'Color', LineColor, 'FontSize', 12)

if params.CBar_on
    pos = get(gca, 'pos');
    colormap(CMAP)
    switch params.view
        case {'lat', 'med'}
            h2 = colorbar(gca, 'Color', LineColor, 'Location', 'southoutside');
            h2.Position = [pos(1)+pos(3)/4, pos(2), pos(3)/2, 0.035];
            set(h2, 'Ticks', params.cbticks, 'TickLabels', params.cblabels);
        case {'post', 'dorsal'}
            h2 = colorbar(gca, 'Color', LineColor, 'Location', 'eastoutside');
            set(h2, 'Ticks', params.cbticks, 'TickLabels', params.cblabels);
    end
    set(h2,'Limits',[0 1]); % some colormaps truncate the limits
end
colormap(CMAP);
hF = params.fig_handle;

%
