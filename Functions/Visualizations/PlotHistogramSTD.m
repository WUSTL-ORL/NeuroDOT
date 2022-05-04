function PlotHistogramSTD(info, params, bthresh)

% PLOTHISTOGRAMSTD A histogram of channel standard deviations.
%
%   PLOTHISTOGRAMSTD(info) plots a histogram of the distribution of
%   channel variances for specified groupings. A red bar will show the
%   threshold value.
%
%   PLOTHISTOGRAMSTD(info, params) allows the user to specify parameters
%   for plot creation.
%
%   PLOTHISTOGRAMSTD(info, params, bthresh) allows the user to specify a
%   threshold value. If no input is supplied, the default value of 0.075
%   will be used.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [200, 200, 560, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       dimension   '2D'                    Dimension of pair radii used.
%       rlimits     (all R2D)               Limits of pair radii displayed.
%       Nnns        (all NNs)               Number of NNs displayed.
%       Nwls        (all WLs)               Number of WLs displayed.
%       xlimits     [0, 50]                 Limits of x-axis.
%       ylimits     [0, 200]                Limits of y-axis.
%       bins        [0:0.5:100]             Number or array of histogram
%                                           bins.
%
% Dependencies: PLOTHISTOGRAMDATA, ISTABLEVAR.
% 
% See Also: FINDGOODMEAS, PLOTCAPGOODMEAS, HISTOGRAM.
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
FieldColor = [0.1, 0.1, 0.1];
cs = unique(info.pairs.WL); % WLs.

pctage_scaling = 100;
Nm = height(info.pairs);
use_NNx_RxD = 'RxD';

lcell = {};
h = [];

if ~exist('bthresh', 'var')  ||  isempty(bthresh)
    bthresh = 0.075;
end
if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [200, 200, 560, 420];
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
if ~isfield(params, 'xlimits')  ||  isempty(params.xlimits)
    params.xlimits = [0, 50]; % STD (as a %)
end
if ~isfield(params, 'ylimits')  ||  isempty(params.ylimits)
    params.ylimits = [0, 200]; % Meas per bin.
end
if ~isfield(params, 'bins')  ||  isempty(params.bins)
    params.bins = [0:0.5:100];
end

if isfield(info.pairs, 'lambda')
    lambdas = unique(info.pairs.lambda);
    lambda_unit1 = '';
    lambda_unit2 = ' nm';
else
    lambdas = cs;
    lambda_unit1 = '\lambda';
    lambda_unit2 = '';
end

%% Plot Data.
hold on
set(gca, 'ColorOrder', [...
    0, 0, 1;... %Blue
    1, 0, 0;... %Red
    1, 1, 1;... %White
    0, 1, 0;... %Green
    1, .5, 0;... %Orange
    1, 0, 1]); % Purple
for k = params.Nwls
    for l = lvar
        switch use_NNx_RxD
            case 'NNx'
                keep = (info.pairs.WL == k)...
                    &  (info.pairs.NN == l);
                % Omitting keepr here because if rlimits present, will
                % never use this branch.
                lcell{end+1} = ['NN', num2str(l), ', ', lambda_unit1,...
                    num2str(lambdas(k)), lambda_unit2];
            case 'RxD'
                keep = (info.pairs.WL == k)...
                    &  (info.pairs.(['r', lower(params.dimension)]) >= params.rlimits(l, 1))...
                    &  (info.pairs.(['r', lower(params.dimension)]) <= params.rlimits(l, 2));
                lcell{end+1} = ['r', lower(params.dimension), ' \in [',...
                    num2str(params.rlimits(l, 1)), ', ',...
                    num2str(params.rlimits(l, 2)), '] mm, ',...
                    lambda_unit1, num2str(lambdas(k)), lambda_unit2];
            case 'all'
                keep = (info.pairs.WL == k);
                lcell{end+1} = ['All r', lower(params.dimension),...
                    ' and NN, ', lambda_unit1, num2str(lambdas(k)), lambda_unit2];
        end
        
        STD = info.MEAS.STD(keep) * pctage_scaling;
        
        if ~isempty(STD)
            h(end+1) = PlotHistogramData(STD, params);
        else
            lcell(end) = [];
        end
    end
end

%% Add threshold.
plot(ones(1,2) * pctage_scaling * bthresh, ylim, 'r', 'LineWidth', 2)

%% Label.
tcell{1} = 'Good Measurements';
if strcmp(use_NNx_RxD, 'all')
    tcell{1} = [tcell{1}, ', All r', lower(params.dimension), ' and NN'];
end
title(tcell, 'Color', LineColor)

ylabel('# of Measurements', 'Color', LineColor)
xlabel('Standard Deviation (%)', 'Color' ,LineColor)

legend(h, lcell, 'Color', FieldColor, 'TextColor', LineColor)



%
