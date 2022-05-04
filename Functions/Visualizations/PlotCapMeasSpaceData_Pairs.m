function PlotCapMeasSpaceData_Pairs(data, info, params)
%
%   PlotCapMeasSpaceData_Pairs(data, info) plots a visualization of a
%   measurements x 1 array of data on a cap layout based on the metadata in
%   info.optodes.  The data value corresponding to each measurement
%   (source-detector pair) is displayed as a line whose color is determined
%   by the data and colormap.  By default, params.mode is set to 'good' so
%   that the plot will only display lines for good measurements in
%   info.MEAS.GI (e.g., determined by FINDGOODMEAS).  If info.MEAS.GI does
%   not exist, then all measurements are treated as good.
%
%   The plot title provides tallies for all specified groupings.  If
%   params.markBadOptodes is set to 1 (the default), the plot title also
%   lists the total number of optodes for which only 33% of their
%   measurements are good, and these optodes are surrounded with white
%   circles.
%
%   PlotCapMeasSpaceData_Pairs(data, info, params) allows the user to
%   specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       dimension   '2D'                    Specifies either a 2D or 3D
%                                           plot rendering.
%       rlimits     (all R2D)               Limits of pair radii displayed.
%       Nnns        (all NNs)               Number of NNs displayed.
%       Nwls        (all WLs)               Number of WLs averaged and
%                                           displayed.
%       mode        'good'                  Display mode. 'good' displays
%                                           channels above noise threhsold,
%                                           'bad' below.
%
% Dependencies: PLOTCAPDATA, ISTABLEVAR.
% 
% See Also: FINDGOODMEAS, PLOTCAP, PLOTCAPMEANLL.
% 
% Copyright (c) 2017 Washington University 
% Authors: Adam T. Eggebrecht, Zachary E. Markow
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

%% Paramters and Initialization.
%BkgdColor = 'k';
LineColor = 'w';
SMarkerColor = [1, .75, .75];
DMarkerColor = [.55, .55, 1];
use_NNx_RxD = 'RxD';

Nm = height(info.pairs);
Ns = numel(unique(info.pairs.Src));
Nd = numel(unique(info.pairs.Det));
cs = unique(info.pairs.WL);
Nc = numel(cs);

thr = 0.33;
N_bad_SD = 0; % Number above threshold.
thr_NNx_RxD = 33;

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params,'title')
    params.title = 'Data on Cap';
end

if (~isfield(params,'markBadOptodes')) || isempty(params.markBadOptodes)
    params.markBadOptodes = 1;
end

if ~isfield(params,'BkgdColor')
    params.BkgdColor = [0,0,0];
end

BkgdColor = params.BkgdColor;

if ~isfield(params,'CBar_on')
    params.CBar_on = 1;
end
if (~isfield(params,'CBar_location')) && params.CBar_on
    params.CBar_location = 'EastOutside';
end

if ~isfield(params, 'mode')  ||  isempty(params.mode)
    params.mode = 'good';
end
if ~isfield(info, 'MEAS')  ||  (isfield(info, 'MEAS')  && ...
        ~isfield(info.MEAS, 'GI'))
    GM = ones(Nm, 1);
else
    GM = info.MEAS.GI;
end
if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D';
end
if (~isfield(params, 'rlimits')  ||  isempty(params.rlimits))  &&...
        (~isfield(params, 'Nnns')  ||  isempty(params.Nnns))
    % If both empty, use ALL.
    use_NNx_RxD = 'all';
    params.rlimits = [min(info.pairs.(['r', lower(params.dimension)])),...
        max(info.pairs.(['r', lower(params.dimension)]))];
    lvar = 1;
else % Otherwise, set defaults if one or the other is missing.
    if ~isfield(params, 'rlimits')  ||  isempty(params.rlimits)
        use_NNx_RxD = 'NNx';
        lvar = params.Nnns;
        thr_NNx_RxD = 2;
    end
    if ~isfield(params, 'Nnns')  ||  isempty(params.Nnns)
        lvar = 1:size(params.rlimits, 1);
    end
end
switch params.dimension
    case '2D'
        Ndim = 2;
        spos = info.optodes.spos2;
        dpos = info.optodes.dpos2;
        MarkerSize = 25;
        if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
            params.fig_size = [20, 200, 1240, 420];
        end
    case '3D'
        Ndim = 3;
        spos = info.optodes.spos3;
        dpos = info.optodes.dpos3;
        MarkerSize = 20;
        if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
            params.fig_size = [20, 200, 560, 560];
        end
end
switch params.mode
    case {'good','both'}
        if strcmp(params.mode,'good')
            modeGM = GM;
        else
            modeGM = true(size(GM));
        end
        %SDLineColor = 'g';
        if (~isfield(params,'line_thickness')) || isempty(params.line_thickness)
            l_Thickness = 0.5 * ones(numel(lvar), 1);
        else
            l_Thickness = repmat(params.line_thickness, numel(lvar), 1);
        end
    case 'bad'
        modeGM = ~GM;
        %SDLineColor = 'r';
        if (~isfield(params,'line_thickness')) || isempty(params.line_thickness)
            l_Thickness = 0.5 + 1.5 * numel(lvar):-1.5:0.5;
        else
            l_Thickness = repmat(params.line_thickness, numel(lvar), 1);
        end
end
% Keep this here, needs fig size from above, and need to spawn figure for
% plotting lines onto.
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, ...
        'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end

%% Plot Src and Det #s.
hold on
params2 = params; params2.mode = 'text';
SrcRGB = repmat(SMarkerColor, Ns, 1); DetRGB = repmat(DMarkerColor, Nd, 1);
tcell = {};
PlotCapData(SrcRGB, DetRGB, info, params2);
switch params.mode
    case {'good','bad'}
        tcell{1} = [upper(params.mode(1)), params.mode(2:end), ' Measurements'];
    case {'both'}
        tcell{1} = 'Good and Bad Measurements';
end
%tcell{1} = params.title;
% switch params.mode
%     case {'good','bad'}
%         PlotCapData(SrcRGB, DetRGB, info, params2);
%         %tcell{1} = [upper(params.mode(1)), params.mode(2:end), ' Measurements'];
%         tcell{end+1} = '';
%         
%     case {'both'}
%         %tcell{1} = ['Both Good and Bad Measurements'];
%         tcell{end+1} = '';
%         params.mode='bad';
%         PlotCapGoodMeas(info, params);
% end

%% Plot GM Lines.

[dataRGB, CMAP, params] = applycmap(data, [], params);

for l = lvar
    % Ignore WLs.
    switch use_NNx_RxD
        case 'RxD'
            keep_NNx_RxD = (info.pairs.(['r', lower(params.dimension)]) >= ...
                params.rlimits(l, 1))...
                &  (info.pairs.(['r', lower(params.dimension)]) <= ...
                params.rlimits(l, 2));
        case 'NNx'
            keep_NNx_RxD = info.pairs.NN == l;
        case 'all'
            keep_NNx_RxD = ones(Nm, 1);
    end
    
    keep = keep_NNx_RxD  &  modeGM;
    keep_idx = find(keep); % NEED THIS FOR THE NEXT SECTION
    
    nMeasKeep = size(keep_idx,1);
    switch params.dimension
        case '2D'
            for keptMeasNum = 1:nMeasKeep
                xCoords = [spos(info.pairs.Src(keep_idx(keptMeasNum)),1); dpos(info.pairs.Det(keep_idx(keptMeasNum)),1)];
                yCoords = [spos(info.pairs.Src(keep_idx(keptMeasNum)),2); dpos(info.pairs.Det(keep_idx(keptMeasNum)),2)];
                plot(xCoords, yCoords, 'Color', dataRGB(keep_idx(keptMeasNum),:), 'LineWidth', l_Thickness(l));
            end
        case '3D'
            for keptMeasNum = 1:nMeasKeep
                xCoords = [spos(info.pairs.Src(keep_idx(keptMeasNum)),1); dpos(info.pairs.Det(keep_idx(keptMeasNum)),1)];
                yCoords = [spos(info.pairs.Src(keep_idx(keptMeasNum)),2); dpos(info.pairs.Det(keep_idx(keptMeasNum)),2)];
                zCoords = [spos(info.pairs.Src(keep_idx(keptMeasNum)),3); dpos(info.pairs.Det(keep_idx(keptMeasNum)),3)];
                plot3(xCoords, yCoords, zCoords, 'Color', dataRGB(keep_idx(keptMeasNum),:), 'LineWidth', l_Thickness(l));
            end
    end
    
    N_GMs(l) = nnz(keep);
    N_Tots(l) = nnz(keep_NNx_RxD);
    
end

%% Calculate and Mark Any Srcs and Dets Above Threshold.

if params.markBadOptodes

    % IE, for each Src or Det, we count it as GM if >1 of its WLx is GM. Since
    % each Src/Det will have WLx * Det # of channels, we can simply count the
    % number of unique Det/Src for each Src/Det's GMs.
    switch use_NNx_RxD
        case {'RxD', 'all'}
            keep_NNx_RxD = info.pairs.(['r', lower(params.dimension)]) <= ...
                thr_NNx_RxD;
        case 'NNx'
            keep_NNx_RxD = info.pairs.NN <= thr_NNx_RxD;
    end
    for s = 1:Ns
        keep = (info.pairs.Src == s)  &  keep_NNx_RxD  &  GM;
        N_GM_Src = numel(unique(info.pairs.Det(keep)));

        % Total number of meas's for this Src.
        tot_Srcs = numel(unique(info.pairs.Det((info.pairs.Src == s)  & ...
            keep_NNx_RxD)));
        if (N_GM_Src / tot_Srcs) < (1 - thr)
            switch params.dimension
                case '2D'
                    plot(spos(s, 1), spos(s, 2), 'ow', 'MarkerSize', MarkerSize)
                case '3D'
                    plot3(spos(s, 1), spos(s, 2), spos(s, 3), 'Marker', 'o',...
                        'MarkerSize', MarkerSize, 'MarkerEdgeColor', LineColor)
            end
            N_bad_SD = N_bad_SD + 1;
        end
    end
    for d = 1:Nd
        keep = (info.pairs.Det == d)  &  keep_NNx_RxD  &  GM;
        N_GM_Det = numel(unique(info.pairs.Src(keep)));

        tot_Dets = numel(unique(info.pairs.Src((info.pairs.Det == d)  & ...
            keep_NNx_RxD)));
        if (N_GM_Det / tot_Dets) < (1 - thr)
            switch params.dimension
                case '2D'
                    plot(dpos(d, 1), dpos(d, 2), 'ow', 'MarkerSize', MarkerSize)
                case '3D'
                    plot3(dpos(d, 1), dpos(d, 2), dpos(d, 3), 'Marker', 'o',...
                        'MarkerSize', MarkerSize, 'MarkerEdgeColor', LineColor)
            end
            N_bad_SD = N_bad_SD + 1;
        end
    end

end

%% Add Title
%tcell{2} = [upper(params.mode(1)), params.mode(2:end), ' Measurements'];
tcell{end+1} = '';
for l = lvar
    tcell{end} = [tcell{end}, num2str(N_GMs(l)), '/', num2str(N_Tots(l)),...
        ' (', num2str(100 * (N_GMs(l) / N_Tots(l)), '%2.0f'), '%) '];
    switch use_NNx_RxD
        case 'NNx'
            tcell{end} = [tcell{end}, 'NN', num2str(l)];
        case 'RxD'
            tcell{end} = [tcell{end}, 'r', lower(params.dimension),...
                ' \in [', num2str(params.rlimits(l, 1)), ', ',...
                num2str(params.rlimits(l, 2)), '] mm'];
    end
    if l ~= lvar(end)
        tcell{end} = [tcell{end}, ', '];
    end
end
if params.markBadOptodes
%     tcell{end+1} = [num2str(N_bad_SD), ' Srcs or Dets have >', ...
%         num2str(100 * thr), '% Bad Measurements'];
    tcell{end+1} = ['>' num2str(100 * thr), '% of meas are bad for ', num2str(N_bad_SD), ' srcs or dets'];
end
% switch use_NNx_RxD
%     case 'NNx'
%         tcell{end} = [tcell{end}, ' \leq NN', num2str(thr_NNx_RxD)];
%     case 'RxD'
%         tcell{end} = [tcell{end}, ', r', lower(params.dimension), ...
%             ' \leq ', num2str(thr_NNx_RxD)];
% end

if ~isempty(params.title)
    % Place user-specified or default title at the top.
    tcell = cat(2, {''}, tcell(1:end));
    tcell{1} = params.title;
end

title(tcell, 'Color', LineColor, 'FontSize', 11)


%% Add colorbar if the user would like.

if params.CBar_on
    
    % Set colormap.
    colormap(CMAP);
    
    % Prepare colorbar labels similarly to PlotSlices.
    params = PrepCustomColorbarLabelParams(params);
    
    % Plot colorbar.
    CB = colorbar('Color', LineColor, 'Location', params.CBar_location);
    set(CB, 'Ticks', params.cbticks);
    set(CB, 'TickLabels', params.cblabels);
    
end


end  % End of function.
