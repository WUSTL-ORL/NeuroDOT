function data=PlotCapMeanLL(data, info, params)

% PLOTCAPMEANLL A visualization of mean light levels overlaid on a cap
% grid.
%
%   PLOTCAPMEANLL(data, info) plots an intensity map of the mean light
%   levels for specified measurement groupings of each optode on the cap
%   and arranges them based on the metadata in "info.optodes".
%
%   PLOTCAPMEANLL(data, info, params) allows the user to specify parameters
%   for plot creation.
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
%       useGM       0                       Use Good Measurements.
%       Cmap.P      'hot'                   Default color mapping.
%
% Dependencies: PLOTCAPDATA, ISTABLEVAR, APPLYCMAP.
%
% See Also: PLOTCAP, PLOTCAPGOODMEAS.
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
LineColor = 'w';
BkgdColor = 'k';
Nm = size(info.pairs,1);
Ns = length(unique(info.pairs.Src));
Nd = length(unique(info.pairs.Det));
cs = unique(info.pairs.WL); % WLs.
use_NNx_RxD = 'RxD';

dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);
Nm = prod(dims(1:end-1));

keep_NNx_RxD = zeros(Nm, 1);

% Hardcoded params for color mapping to the units scaling used here.
scale_order = 1e-2;
base = 1e9;
dr = 3;
params.mode = 'patch';
if ~isfield(params,'Scale'), params.Scale=dr;end
if ~isfield(params,'abs'), params.abs=0;end % Absolute or relative coloring
params.PD = 1;
params.Th.P = 0;
params.DR = 1000;

if (~isfield(params, 'Cmap')  ||  isempty(params.Cmap))...
        ||  (~isfield(params.Cmap, 'P')  ||  isempty(params.Cmap.P))
    params.Cmap.P = 'hot';
end
if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D';
end
if (~isfield(params, 'rlimits')  ||  isempty(params.rlimits))  &&...
        (~isfield(params, 'Nnns')  ||  isempty(params.Nnns))
    % If both empty, use ALL.
    use_NNx_RxD = 'all';
    lvar = 1;
else % Otherwise, set defaults if one or the other is missing.
    if ~isfield(params, 'rlimits')  ||  isempty(params.rlimits)
        use_NNx_RxD = 'NNx';
        lvar = params.Nnns;
    end
    if ~isfield(params, 'Nnns')  ||  isempty(params.Nnns)
        lvar = 1:size(params.rlimits, 1);
    end
end
if ~isfield(params, 'Nwls')  ||  isempty(params.Nwls)
    params.Nwls = cs';
end
if ~isfield(params, 'useGM')  ||  isempty(params.useGM)
    params.useGM = 0;
end
if ~params.useGM  ||  ~isfield(info, 'MEAS')  ||  (isfield(info, 'MEAS')...
        &&  ~isfield(info.MEAS, 'GI'))
    GM = ones(Nm, 1);
else
    GM = info.MEAS.GI;
end
if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D'; % '2D' | '3D'
end
if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    switch params.dimension
        case '2D'
            params.fig_size = [20, 200, 1240, 420];
        case '3D'
            params.fig_size = [20, 200, 560, 560];
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

if isfield(info.pairs, 'lambda')
    lambdas = unique(info.pairs.lambda);
    if numel(params.Nwls) > 1
        lambda_unit1 = '[';
        lambda_unit2 = '] nm';
    else
        lambda_unit1 = '';
        lambda_unit2 = ' nm';
    end
else
    lambdas = cs;
    lambda_unit1 = '\lambda';
    lambda_unit2 = '';
end

%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Average LLs of each Src and Det.
% Average LLs over time first.
data = mean(data, 2);

switch use_NNx_RxD
    case 'RxD'
        for k = 1:size(params.rlimits, 1)
            keep_NNx_RxD = keep_NNx_RxD...
                |  ((info.pairs.(['r', lower(params.dimension)]) >= ...
                params.rlimits(k, 1))...
                &  (info.pairs.(['r', lower(params.dimension)]) <= ...
                params.rlimits(k, 2)));
        end
    case 'NNx'
        for k = params.Nnns
            keep_NNx_RxD = keep_NNx_RxD  |  (info.pairs.NN == k);
        end
    case 'all'
        keep_NNx_RxD = ~keep_NNx_RxD;
end

% Desired WLs.
keepWL = zeros(Nm, 1);
for k = params.Nwls
    keepWL = keepWL  |  (info.pairs.WL == k);
end

if params.abs
    % Src averages. (in desired R)
    for s = 1:Ns
        keep = (info.pairs.Src == s)  &  keep_NNx_RxD  &  keepWL  &  GM;
        Srcval(s) = log10(base * mean(data(keep)));
    end
    
    % Det averages. (in desired R)
    for d = 1:Nd
        keep = (info.pairs.Det == d)  &  keep_NNx_RxD  &  keepWL  &  GM;
        Detval(d) = log10(base * mean(data(keep)));
    end
    
else
    
    % Src averages. (in desired R)
    for s = 1:Ns
        keep = (info.pairs.Src == s)  &  keep_NNx_RxD  &  keepWL  &  GM;
        Srcval(s) = log10(mean(data(keep)));
    end
    
    % Det averages. (in desired R)
    for d = 1:Nd
        keep = (info.pairs.Det == d)  &  keep_NNx_RxD  &  keepWL  &  GM;
        Detval(d) = log10(mean(data(keep)));
    end
    
end

%% Scaling and Color Mapping
if params.abs % Absolute data values color mapping
Srcval = Srcval - (log10(base) - dr);
Detval = Detval - (log10(base) - dr);
else        % Relative data values
    M=max([max(Srcval(:)),max(Detval)]);
    Srcval = Srcval - (M - dr);
    Detval = Detval - (M - dr);    
end

[SDRGB, CMAP] = applycmap([Srcval(:); Detval(:)], [], params);
SrcRGB = SDRGB(1:Ns, :);
DetRGB = SDRGB(Ns+1:end, :);

%% Send Light Levels to PlotCapData.
PlotCapData(SrcRGB, DetRGB, info, params)

%% Add Colorbar.
colormap(CMAP)
pos = get(gca, 'pos');

if params.abs % Absolute data values color mapping
CB = colorbar('YLim', [0, 1], 'YTick', [0, 0.5, 1],...
    'YTickLabel',{'10nW', '', '1\muW'}, 'Color', LineColor,...
    'Location', 'southoutside');
else
CB = colorbar('YLim', [0, 1], 'YTick', [0, 0.5, 1],...
    'YTickLabel',{[num2str(10^-(dr-2)),'% max'], '', 'max'},...
    'Color', LineColor,'Location', 'southoutside');
    
end

switch params.dimension
    case '2D'
        CB.Position = [pos(1)+pos(3)/3, pos(2), pos(3)/3, pos(4)/35];
    case '3D'
        CB.Position = [pos(1)+pos(3)/3, pos(2)+0.11, pos(3)/3, 0.015];
end

%% Add Title.
tcell{1} = 'Mean Light Levels';
switch use_NNx_RxD
    case 'all'
        tcell{end+1} = ['All r', lower(params.dimension), ' and NN'];
    case 'RxD'
        tcell{end+1} = ['r', lower(params.dimension), ' \in '];
        for l = lvar
            if (l == lvar(1))  &&  (numel(lvar) > 1)
                tcell{end} = [tcell{end}, '('];
            end
            tcell{end} = [tcell{end}, '[',...
                num2str(params.rlimits(l, 1)), ', ',...
                num2str(params.rlimits(l, 2)), ']'];
            if l ~= lvar(end)
                tcell{end} = [tcell{end}, ', '];
            elseif numel(lvar) > 1
                tcell{end} = [tcell{end}, ') mm'];
            else
                tcell{end} = [tcell{end}, ' mm'];
            end
        end
    case 'NNx'
        tcell{end+1} = ['NN', num2str(params.Nnns, '%d')];
end
if isfield(info.pairs, 'lambda')
    if numel(params.Nwls) > 1
        wlstr = num2str(lambdas', '%d, ');
        wlstr(end)=[];
    else
        wlstr = num2str(lambdas(params.Nwls));
    end
else
    wlstr = num2str(cs', '%d');
end
tcell{end} = [tcell{end}, ', ', lambda_unit1, wlstr, lambda_unit2];
if params.useGM
    tcell{end} = [tcell{end}, ', GM'];
end
title(tcell, 'Color', LineColor)



%
