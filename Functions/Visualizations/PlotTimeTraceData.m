function h = PlotTimeTraceData(data, time, params)

% PLOTTIMETRACEDATA A basic time traces plotting function.
%
%   PLOTTIMETRACEDATA(data, time) takes a light-level array "data" of the
%   MEAS x TIME format, and plots its time traces.
%
%   PLOTTIMETRACEDATA(data, time, params) allows the user to specify
%   parameters for plot creation.
%
%   h = PLOTTIMETRACEDATA(...) passes the handles of the plot line objects
%   created.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [200, 200, 560, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       xlimits     'auto'                  Limits of x-axis.
%       xscale      'linear'                Scaling of x-axis.
%       ylimits     'auto'                  Limits of y-axis.
%       yscale      'linear'                Scaling of y-axis.
%
%   Shaded event periods (optional; drawn behind the traces at full plot
%   height). Turn on with event_shade = 1, then give EITHER explicit times
%   OR a paradigm structure:
%       event_shade   0                     Flag: 1 = draw the gray event
%                                           boxes, 0 = don't (the event
%                                           fields below are then ignored).
%       event_onsets  []                    Event onset times, in the same
%                                           units as "time" (e.g., seconds).
%       event_dur     []                    Event duration(s), same units:
%                                           a scalar (e.g., 10 for a 10 s
%                                           block) or one value per onset.
%       paradigm      []                    An info.paradigm structure.
%                                           Onsets are time(synchpts(Pulse_N))
%                                           with N = params.event_pulse.
%                                           Without event_dur, each event
%                                           lasts until the next synchpt.
%                                           "time" must be the full run's
%                                           time vector (synchpts index it).
%       event_pulse   2                     Which Pulse_N to shade.
%       event_color   [0.3, 0.3, 0.3]       Box color (RGB).
%
%   Examples:
%       % 10 s tensing blocks at known onsets (s)
%       params.event_shade = 1;
%       params.event_onsets = tense_onsets; params.event_dur = 10;
%       PlotTimeTraceData(data, t, params);
%
%       % Blocks straight from a NeuroDOT paradigm
%       t = (0:size(data, 2) - 1) / info.system.framerate;
%       params.event_shade = 1; params.paradigm = info.paradigm;
%       PlotTimeTraceData(data, t, params);
%
% See Also: PLOTTIMETRACEALLMEAS, PLOTTIMETRACEMEAN.
%
% Copyright (c) 2017 Washington University
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable,
% non-exclusive, royalty-free, non-commercial, research license to use
% and copy the computer code that is provided here (the Software).
% You agree to include this license and the above copyright notice in
% all copies of the Software.  The Software may not be distributed,
% shared, or transferred to any third party.  This license does not
% grant any rights or licenses to any other patents, copyrights, or
% other forms of intellectual property owned or controlled by Washington
% University.
%
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR
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

dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);

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
if ~isfield(params, 'xlimits')  ||  isempty(params.xlimits)
    params.xlimits = 'auto';
end
if ~isfield(params, 'xscale')  ||  isempty(params.xscale)
    params.xscale = 'linear';
end
if ~isfield(params, 'ylimits')  ||  isempty(params.ylimits)...
        ||  (all(data(:) > params.ylimits(2))...
        &&  all(data(:) < params.ylimits(1)))   % data(:): also works for N-D input
    params.ylimits = 'auto';
end
if ~isfield(params, 'yscale')  ||  isempty(params.yscale)
    params.yscale = 'linear';
end
if ~isfield(params, 'event_shade')  ||  isempty(params.event_shade)
    params.event_shade = 0;
end
if ~isfield(params, 'event_onsets')
    params.event_onsets = [];
end
if ~isfield(params, 'event_dur')
    params.event_dur = [];
end
if ~isfield(params, 'paradigm')
    params.paradigm = [];
end
if ~isfield(params, 'event_pulse')  ||  isempty(params.event_pulse)
    params.event_pulse = 2;
end
if ~isfield(params, 'event_color')  ||  isempty(params.event_color)
    params.event_color = [0.3, 0.3, 0.3];
end

if ~exist('time', 'var')
    time = 1:Nt;
end

%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Plot Data.
if isempty(data),data=zeros(size(time));end
h = plot(time, data', 'LineWidth', 1); % Plot mean of measurements vs. time.

%% Format Plot.
box on

a = gca;
a.Color = BkgdColor;
a.XColor = LineColor;
a.YColor = LineColor;

%% Scaling.
ylim(params.ylimits)
xlim(params.xlimits)
a.YScale = params.yscale;
a.XScale = params.xscale;

%% Shaded event periods.
onsets = params.event_onsets(:);
durs = params.event_dur(:);
if ~params.event_shade           % flag off: no boxes
    onsets = [];
elseif isempty(onsets)  &&  ~isempty(params.paradigm)
    % Onsets (and, if needed, durations) from a NeuroDOT paradigm structure
    sp = params.paradigm.synchpts(:);
    idx = params.paradigm.(['Pulse_', num2str(params.event_pulse)])(:);
    if any(sp(idx) > numel(time))
        error(['*** paradigm synchpts exceed the length of "time". Pass the ', ...
            'full run, or use params.event_onsets instead. ***'])
    end
    onsets = reshape(time(sp(idx)), [], 1);
    if isempty(durs)
        nxt = min(idx + 1, numel(sp));
        ends = reshape(time(min(sp(nxt), numel(time))), [], 1);
        ends(idx == numel(sp)) = time(end);
        durs = ends - onsets;
    end
end
if ~isempty(onsets)
    if isempty(durs)
        error('*** params.event_dur is required with params.event_onsets. ***')
    end
    if isscalar(durs)
        durs = repmat(durs, size(onsets));
    elseif numel(durs) ~= numel(onsets)
        error('*** params.event_dur must be a scalar or match the number of onsets. ***')
    end
    yl = a.YLim;                 % current limits (also when 'auto')
    held = ishold(a);
    hold(a, 'on')
    for k = 1:numel(onsets)
        p = patch(a, [onsets(k), onsets(k) + durs(k), onsets(k) + durs(k), onsets(k)], ...
            [yl(1), yl(1), yl(2), yl(2)], params.event_color, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(p, 'bottom');    % traces stay on top of the boxes
    end
    ylim(a, yl)                  % boxes must not change the y-limits
    if ~held, hold(a, 'off'); end
end



%
