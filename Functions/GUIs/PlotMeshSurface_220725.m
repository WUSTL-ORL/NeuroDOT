function PlotMeshSurface_200507(mesh, params)

% PLOTMESHSURFACE Creates a 3D surface mesh visualization.
% 
%   PLOTMESHSURFACE(mesh) creates a 3D visualization of the surface mesh
%   "mesh". If no region data is provided in "mesh.region", all nodes will
%   be assumed to form a single region. If a field "data" is provided as
%   part of the "mesh" structure, that data will be used to color the
%   visualization. If both "data" and "region" are present, the "region"
%   values are used as an underlay for the colormapping.
%   
% 
%   PLOTMESHSURFACE(mesh, params) allows the user to specify parameters
%    for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]        Default figure position
%                                               vector.
%       fig_handle  (none)                      Specifies a figure to
%                                               target.
%       Cmap.P      'jet'                       Colormap for positive data
%                                               values.
%       BG          [0.8, 0.8, 0.8]             Background color, as an RGB
%                                               triplet.
%       orientation 't'                         Select orientation of
%                                               volume. 't' for transverse,
%                                               's' for sagittal, 'coord'
%                                               for simple coordinates
%                                               (default) with no flips
%       cboff       0                           If set to 1, no colorbar is
%                                               displayed
% 
%   Note: APPLYCMAP has further options for using "params" to specify
%   parameters for the fusion, scaling, and colormapping process.
% 
% Dependencies: APPLYCMAP
% 
% See Also: PLOTSLICES, PLOTCAP, CAP_FITTER.
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


%% Parameters and Initialization
LineColor = 'k';

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'orientation')
    params.orientation = 'coord'; % need this
end

if ~isfield(params, 'BGC')  ||  isempty(params.BGC)
    params.BGC = [1,1,1]; % Background color of figure, need this
end

% if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size) % don't need this
%     params.fig_size = [20, 200, 560, 560];
% end
% if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle) % dont need new figure
%     params.fig_handle = figure('Color',  'w',...
%         'Position', params.fig_size);
%     new_fig = 1;
% else
%     switch params.fig_handle.Type
%         case 'figure'
%             set(groot, 'CurrentFigure', params.fig_handle);
%         case 'axes'
%             set(gcf, 'CurrentAxes', params.fig_handle);
%     end
% end

fig_handle = params.fig_handle;

if ~isfield(params,'Cmap'), params.Cmap=struct;end
if ~isstruct(params.Cmap)
    temp = params.Cmap;
    params.Cmap = [];
    params.Cmap.P = temp; % params.Cmap.P = 'jet'
end
if ~isfield(params.Cmap,'P'), params.Cmap.P='gray';end
if ~isfield(params,'FaceAlpha'), params.FaceAlpha=1;end % Transparency
if ~isfield(params,'EdgeAlpha'), params.EdgeAlpha=1;end % Transparency
if ~isfield(params,'OL'), params.OL=0;end
if ~isfield(params,'reg'), params.reg=1;end
if ~isfield(params,'TC'),params.TC=0;end  
if ~isfield(params,'PD'), params.PD=0;end
if ~isfield(params,'FaceColor'), params.FaceColor=[0.25, 0.25, 0.25];end
if ~isfield(params,'EdgeColor'), params.EdgeColor='k';end
if ~isfield(params,'EdgesON'), params.EdgesON=1;end
if isfield(params,'AmbientStrength'), 
    AmbientStrength=params.AmbientStrength;
else
    AmbientStrength=0.25; % need this
end
if isfield(params,'DiffuseStrength'), 
    DiffuseStrength=params.DiffuseStrength;
else
    DiffuseStrength=0.75; % need this
end
if isfield(params,'SpecularStrength'), 
    SpecularStrength=params.SpecularStrength;
else
    SpecularStrength=0.1; % need this
end
if ~params.EdgesON, params.EdgeColor='none';end
if ~isfield(params, 'cbmode')  ||  isempty(params.cbmode)
    params.cbmode = 0; % need this
end
if ~isfield(params, 'cboff')  
    params.cboff = 0;
end

if ~isfield(params, 'ctx')  ||  isempty(params.ctx)
    params.ctx = 'std';
end
switch params.ctx
    case 'std'
        % do nothing , need this
    case 'inf'
        mesh.nodes=mesh.Inodes;
    case 'vinf'
        mesh.nodes=mesh.VInodes;
    case 'MTH'
        mesh.nodes=mesh.MTHnodes;
    case 'SPH'
        mesh.nodes=mesh.SPHnodes;
end



%% Get face centers of elements for S/D pairs.
switch size(mesh.elements, 2)
    case 4  % extract surface mesh from volume mesh, need this
        TR = triangulation(mesh.elements, mesh.nodes);
        [m.elements, m.nodes] = freeBoundary(TR);
        [~, Ib] = ismember(m.nodes, mesh.nodes, 'rows');
        Ib(Ib == 0) = []; % Clear zero indices.
        if isfield(mesh,'region'), m.region=mesh.region(Ib);end % don't need 
        if isfield(mesh,'data'), m.data=mesh.data(Ib);end % don't need
    case 3
        m=mesh;
end


%% mesh.data and mesh.region together determine coloring rules

if isfield(m,'region') && ~params.reg % don't need
    m=rmfield(m,'region');
end
    
if ~isfield(m,'data')       % NO DATA, need this
    if ~isfield(m,'region') % no data, no regions, need this
        cb=0;
        FaceColor = params.FaceColor;
        EdgeColor = params.EdgeColor;
        FaceLighting = 'flat';
%         AmbientStrength = 0.5;        
        h = patch(fig_handle,'Faces', m.elements, 'Vertices',m.nodes,...
            'EdgeColor', EdgeColor,'EdgeAlpha',params.EdgeAlpha,... 
            'FaceColor', FaceColor,...
            'FaceLighting', FaceLighting,'FaceAlpha',params.FaceAlpha,...
            'AmbientStrength', AmbientStrength,'DiffuseStrength',... 
            DiffuseStrength,'SpecularStrength', SpecularStrength);  %% this section draws patches on figure. convert patches to fit on UI Axes      
        
    else                      % data are regions               
        params.PD=1;
        params.TC=1;
        if max(m.region(:))>1
            params.DR=max(m.region(:));
            tempCmap=params.Cmap.P;
            if ischar(params.Cmap.P)
            params.Cmap.P=eval([tempCmap, '(', num2str(params.DR), ');']);
            end
        else
            tempCmap=params.Cmap.P;
            if ischar(params.Cmap.P)
            params.Cmap.P=eval([tempCmap, '(', num2str(2), ');']);
            end
            if size(params.Cmap.P,1)>1
            params.Cmap.P=max(params.Cmap.P,1);
            end
            params.DR=1;
        end
        cb=1;
        CMAP=params.Cmap.P;
        EdgeColor = params.EdgeColor;
        FaceColor = 'flat';
        FaceLighting = 'gouraud';
%         AmbientStrength = 0.25;
%         DiffuseStrength = 0.75; % or 0.75
%         SpecularStrength = 0.1;        
        FV_CData=params.Cmap.P(mode(reshape(m.region(m.elements(:)),[],3),2),:);        
        h = patch(fig_handle,'Faces', m.elements, 'Vertices', m.nodes,...
            'EdgeColor', EdgeColor, 'FaceColor', FaceColor,...
            'FaceVertexCData', FV_CData, 'FaceLighting', FaceLighting,...
            'AmbientStrength', AmbientStrength, 'DiffuseStrength',... 
            DiffuseStrength,'SpecularStrength', SpecularStrength);        
    end
    
else                        % DATA
    if ~isfield(m,'region') % no regions
        if ~isfield(params,'BG'),params.BG=[0.25, 0.25, 0.25];end
        cb=2;
        [FV_CData,CMAP] = applycmap(m.data, [], params);
    else                    % with regions: grayscale underlay
        cb=3;
        [FV_CData,CMAP] = applycmap(m.data, m.region, params);
    end
    EdgeColor = params.EdgeColor;
    FaceColor = 'interp';
    FaceLighting = 'gouraud';
%     AmbientStrength = 0.25;
%     DiffuseStrength = 0.75; % or 0.75
%     SpecularStrength = 0.5;
    h = patch(fig_handle,'Faces', m.elements, 'Vertices', m.nodes,...
        'EdgeColor', EdgeColor, 'FaceColor', FaceColor,...
        'FaceVertexCData', FV_CData, 'FaceLighting', FaceLighting,...
        'AmbientStrength', AmbientStrength, 'DiffuseStrength',...
        DiffuseStrength,'SpecularStrength', SpecularStrength);
end
        
        
set(fig_handle, 'Color', params.BGC);%, 'XTick', [], 'YTick', [], 'ZTick', []);

switch params.orientation
    case 's'
        set(fig_handle, 'ZDir', 'rev');
    case 't'
        set(fig_handle, 'XDir', 'rev');
    case 'c'
        set(fig_handle, 'YDir', 'rev');
end

axis(fig_handle,'image')
% axis off
hold(fig_handle,'on')
rotate3d(fig_handle,'on')


%% Set additional lighting
% warning('off','all')
% % Lower lighting
% light(fig_handle,'Position', [-140, 90, -100], 'Style', 'local')
% light(fig_handle,'Position', [-140, -350, -100], 'Style', 'local')
% light(fig_handle,'Position', [300, 90, -100], 'Style', 'local')
% light(fig_handle,'Position', [300, -350, -100], 'Style', 'local')
% 
% % Higher lighting
% light(fig_handle,'Position', [-140, 90, 360], 'Style', 'local');
% light(fig_handle,'Position', [-140, -350, 360], 'Style', 'local');
% light(fig_handle,'Position', [300, 90, 360], 'Style', 'local');
% light(fig_handle,'Position', [300, -350, 360], 'Style', 'local');

% xlabel(fig_handle,'X', 'Color', LineColor)  % changing x labels to black
% ylabel(fig_handle,'Y', 'Color', LineColor)
% zlabel(fig_handle,'Z', 'Color', LineColor)

% if new_fig
% %     view(163, -86)
% % end
% 
% if isfield(params,'side') % don't need
%     switch params.side
%         case 'post'
%             light('Position',[-500,-20,0],'Style','local');
%             light('Position',[500,-20,0],'Style','local');
%             light('Position',[0,-200,50],'Style','local');
%         case 'dorsal'
%             light('Position',[-500,-20,100],'Style','local');
%             light('Position',[500,-20,100],'Style','local');
%             light('Position',[100,-200,200],'Style','local');
%             light('Position',[100,200,200],'Style','local');
%             light('Position',[0,200,0],'Style','local');
%             light('Position',[200,200,0],'Style','local');
%             light('Position',[100,500,100],'Style','local');
%         case 'coronal'
%             if ~any(m.nodes(:)<0)
%                 light('Position',[-500,-20,100],'Style','local');
%                 light('Position',[500,-20,100],'Style','local');
%                 light('Position',[100,-200,200],'Style','local');
%                 light('Position',[100,200,200],'Style','local');
%                 light('Position',[0,200,0],'Style','local');
%                 light('Position',[200,200,0],'Style','local');
%                 light('Position',[100,500,100],'Style','local');
%             else
%                 mm=mean(m.nodes);
%                 x=-mm(1);y=mm(2);z=-mm(3);
%                 light('Position',[-(x+400),y-150,z],'Style','local');
%                 light('Position',[x+400,y-150,z],'Style','local');
%                 light('Position',[x,-(y+50),z+100],'Style','local');
%                 light('Position',[x,y+50,z+100],'Style','local');
%                 % light('Position',[100,200,0],'Style','local');
%                 light('Position',[x-100,y+50,z-100],'Style','local');
%                 light('Position',[x+100,y+50,z-100],'Style','local');
%                 light('Position',[x,y+350,z],'Style','local');
%             end
%             
%     end
% end


%% Set colorbar properties and add colorbar if desired
if ~params.cboff
if cb
    params.cbmode=1; % cbmode was set to off so don't need
    if isfield(params, 'Scale')  &&  ~isempty(params.Scale)
        c_max = params.Scale;
    else
        switch cb
            case 1          % data are regions
                c_max =params.DR;
                
            case 2          % data with no regions
                c_max =0.8*max(m.data(:));
                
            case 3          % data with regions
                c_max =0.8*max(m.data(:));   
        end        
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
    
    % Add the colorbar.
    colormap(gca,CMAP)
    if ~params.cboff
        h2 = colorbar('Color', LineColor);
        if params.cbmode
            set(h2, 'Ticks', params.cbticks, 'TickLabels', params.cblabels);
        end
    end
end
end


%
