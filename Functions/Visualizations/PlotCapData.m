function PlotCapData(SrcRGB, DetRGB, info, params)

% PLOTCAPDATA A basic plotting function for generating and labeling cap grids.
%
%   PLOTCAPDATA(SrcRGB, DetRGB, info) plots the input RGB information in
%   one of three modes:
% 
%   'text' - Optode numbers are arranged in a cap grid and colored with the
%   RGB input.
%   'patch' - Optodes are plotted as patches and colored with the RGB
%   input.
%   'textpatch' - Optodes are plotted as patches and colored with the RGB
%   input, with optode numbers overlain in white.
%
%   PLOTCAPDATA(SrcRGB, DetRGB, info, params) allows the user to specify
%   parameters for plot creation.
% 
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       dimension   '2D'                    Specifies either a 2D or 3D
%                                           plot rendering.
%       mode        'textpatch'             Display mode.
% 
% See Also: PLOTCAP, PLOTCAPGOODMEAS, PLOTCAPMEANLL.
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

%% Parameters and Initialization.

% Nm = height(info.pairs);
Ns = length(info.optodes.spos3);
Nd = length(info.optodes.dpos3);

Srcs = [1:Ns];%unique(info.pairs.Src);
Dets = [1:Nd];%unique(info.pairs.Det);

BkgdColor = [0, 0, 0]; % KEEP THESE LIKE THIS!
% LineColor = [1, 1, 1]; % They are read as RGB triplets by later functions.

box_buffer = max(info.pairs.r2d)/10;
new_fig = 0;

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D'; % '2D' | '3D'
end

if ~isfield(params,'rhombus'),params.rhombus=1;end
if ~isfield(params,'eeg_style'),params.eeg_style=0;end

if ~isfield(params,'LineColor'),params.LineColor=[1,1,1];end
switch params.dimension
    case '2D'
        % Get optode positions.
        spos = info.optodes.spos2;
        dpos = info.optodes.dpos2;
        l = min(pdist2(info.optodes.spos2,info.optodes.dpos2),[],'all')/2;
        if params.rhombus
            % Calculate side length and create square vectors.
            xsq = [l, 0, -l, 0]; % Okay, they're really rhombi, but you get the point.
            ysq = [0, -l, 0, l];
        else
            % Calculate side length and create square vectors.
%             nn1S=info.pairs.Src(find(info.pairs.NN==1,1));
%             nn1D=info.pairs.Det(find(info.pairs.NN==1,1));
            xsq = [l, -l, -l, l]; % Okay, they're really rhombi, but you get the point.
            ysq = [l, l, -l, -l];
        end
        
        if params.eeg_style
            [xsq,ysq,z]=cylinder(2);
            xsq=xsq(1,:).*(l/2);
            ysq=ysq(1,:).*(l/2);
        end
            
        
        % Default figure size.
        if ~isfield(params, 'fig_handle')
        if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
            params.fig_size = [20, 200, 1240, 420];
        end
        end
    case '3D'
        % Get optode positions.
        spos = info.optodes.spos3;
        dpos = info.optodes.dpos3;
        
%         FieldColor = [0.25, 0.25, 0.25];
%         [xdir, ydir, zdir] = CheckOrientation(info);
        
        % Default figure size.
        if ~isfield(params, 'fig_handle')
        if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
            params.fig_size = [20, 200, 560, 560];
        end
        end
end

if ~isfield(params, 'mode')  ||  isempty(params.mode)
    params.mode = 'textpatch'; % 'textpatch' | 'text' | 'patch'
end
switch params.mode
    case 'text'
        STextColor = SrcRGB;
        DTextColor = DetRGB;
        switch params.dimension
            case '2D'
                TextSize = 10;
                if params.eeg_style == 1
                    STextColor(:, :) = STextColor(:, :)*0.6;
                    DTextColor(:, :)= DTextColor(:, :)*0.6;
                    TextSize = TextSize*0.8;
                end
            case '3D'
                TextSize = 8;
        end
    case 'patch'
        switch params.dimension
            case '2D'
                % Nothing for now.
            case '3D'
                MarkerSize = 6;
                SMarkerEdgeColor = SrcRGB;
                DMarkerEdgeColor = DetRGB;
        end
    case 'textpatch'
        % Smaller text and larger markers needed when plotting text over
        % markers.
        switch params.dimension
            case '2D'
                TextSize = 8;
                STextColor = repmat(params.LineColor, Ns, 1);
                DTextColor = repmat(params.LineColor, Nd, 1);
                
            case '3D'
                TextSize = 6;
                MarkerSize = 9;
                STextColor = repmat(params.LineColor, Ns, 1);
                DTextColor = repmat(params.LineColor, Nd, 1);
                SMarkerEdgeColor = STextColor;
                DMarkerEdgeColor = DTextColor;
        end
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

%% Plot an invisible box around grid extrema.
% This is a hack to ensure proper boundaries are kept. MATLAB does not
% count text objects when auto-scaling axes to data. Thus, when other
% objects are only drawn in specific regions of cap grid, axes only get
% scaled to those objects, and plot looks zoomed in. By adding an invisible
% box around the entire cap grid, MATLAB auto-scales axes to this box.
hold on
switch params.dimension
    case '2D'
        mins = min([spos; dpos]) - box_buffer; maxes = max([spos; dpos]) + box_buffer;
        right = maxes(1); left = mins(1); top = maxes(2); bot = mins(2);
        plot([left; right; right; left; left], [top; top; bot; bot; top],...
            'LineStyle', 'none')
    case '3D'
        mins = min([spos; dpos]) - box_buffer; maxes = max([spos; dpos]) + box_buffer;
        right = maxes(1); left = mins(1); front = maxes(2);
        back = mins(2); top = maxes(3); bot = mins(3);
        plot3([left; right; right; left; left; right; right; left; left],...
            [front; front; front; front; back; back; back; back; front],...
            [top; top; bot; bot; bot; bot; top; top; top], 'LineStyle', 'none')
end

%% (3D) Allow Rotation and Force Aspect Ratio.
if strcmp(params.dimension, '3D')
    view(55, 25)
    rotate3d on
    pbaspect manual
    set(gca,'CameraViewAngle', 16, 'XLimMode', 'manual',...
         'YLimMode', 'manual','ZLimMode', 'manual');        
end

%% Draw data.
switch params.dimension
    case '2D'
        
        if params.eeg_style
            % Reshape input before drawing - this is crucial for patch.
            SrcRGB = reshape(SrcRGB, [], 1, 3);
            DetRGB = reshape(DetRGB, [], 1, 3);
            
            
            % Srcs
            patch([repmat(spos(:, 1), 1, size(xsq,2)) + repmat(xsq, Ns, 1)]',...
                [repmat(spos(:, 2), 1, size(xsq,2)) + repmat(ysq, Ns, 1)]',...
                SrcRGB, 'EdgeColor', 'none')
            % Dets
            patch([repmat(dpos(:, 1), 1, size(xsq,2)) + repmat(xsq, Nd, 1)]',...
                [repmat(dpos(:, 2), 1, size(xsq,2)) + repmat(ysq, Nd, 1)]',...
                DetRGB, 'EdgeColor', 'none')
            
        elseif ~isempty(strfind(params.mode, 'patch'))
            % Reshape input before drawing - this is crucial for patch.
            SrcRGB = reshape(SrcRGB, [], 1, 3);
            DetRGB = reshape(DetRGB, [], 1, 3);
            
            % Srcs
            patch([repmat(spos(:, 1), 1, 4) + repmat(xsq, Ns, 1)]',...
                [repmat(spos(:, 2), 1, 4) + repmat(ysq, Ns, 1)]',...
                SrcRGB, 'EdgeColor', 'none')
            % Dets
            patch([repmat(dpos(:, 1), 1, 4) + repmat(xsq, Nd, 1)]',...
                [repmat(dpos(:, 2), 1, 4) + repmat(ysq, Nd, 1)]',...
                DetRGB, 'EdgeColor', 'none')
        end
        if ~isempty(strfind(params.mode, 'text'))
            for s = 1:Ns
                text(spos(s, 1), spos(s, 2), num2str(Srcs(s)),...
                    'Color', STextColor(s, :), 'HorizontalAlignment', 'center',...
                    'FontSize', TextSize, 'FontWeight', 'bold')
            end
            for d = 1:Nd
                text(dpos(d, 1), dpos(d, 2), num2str(Dets(d)),...
                    'Color', DTextColor(d, :), 'HorizontalAlignment', 'center',...
                    'FontSize', TextSize, 'FontWeight', 'bold')
            end
        end
        axis off
    case '3D'
        if ~isempty(strfind(params.mode, 'patch'))
            % Srcs
            for s = 1:Ns
                plot3(spos(s, 1), spos(s, 2), spos(s, 3), 'Color', SrcRGB(s, :),...
                    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize,...
                    'MarkerFaceColor', SrcRGB(s, :),...
                    'MarkerEdgeColor', SMarkerEdgeColor(s, :));
            end
            % Dets
            for d = 1:Nd
                plot3(dpos(d, 1), dpos(d, 2), dpos(d, 3), 'Color', DetRGB(d, :),...
                    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize,...
                    'MarkerFaceColor', DetRGB(d, :),...
                    'MarkerEdgeColor', DMarkerEdgeColor(d, :));
            end
        end
        if ~isempty(strfind(params.mode, 'text'))
            for s = 1:Ns
                text(spos(s, 1), spos(s, 2), spos(s, 3), num2str(Srcs(s)),...
                    'Color', STextColor(s, :), 'HorizontalAlignment', 'center',...
                    'FontSize', TextSize, 'FontWeight', 'bold')
            end
            for d = 1:Nd
                text(dpos(d, 1), dpos(d, 2), dpos(d, 3), num2str(Dets(d)),...
                    'Color', DTextColor(d, :), 'HorizontalAlignment', 'center',...
                    'FontSize', TextSize, 'FontWeight', 'bold')
            end
        end
end

%% Remove Axes and Resize.
axis image



%
