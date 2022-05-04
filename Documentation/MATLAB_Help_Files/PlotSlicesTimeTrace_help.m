%% PlotSlicesTimeTrace
% Creates an interactive 4D plot.
%
%% Description
% |PlotSlicesTimeTrace(underlay, infoVol, params, overlay)| uses the same
% basic inputs as the function |PlotSlices| to create an interactive 3-axis
% plot of a 4D volume, plus an axis for the time trace of the selected
% voxel.
%
% |PlotSlicesTimeTrace(..., info)| allows the user to input an |info|
% structure to provide information about the 4D volume's native framerate.
% The default framerate is |1 Hz|.
% 
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>fig_size</td><td>[20, 200, 840, 720]</td><td>Default figure
% position vector.</td></tr>
% <tr><td>fig_handle</td><td>(none)</td><td>Specifies a figure to target.
% If empty, spawns a new figure.</td></tr>
% <tr><td>CH</td><td>1</td><td>Turns crosshairs on (1) and off
% (0).</td></tr>
% <tr><td>Scale</td><td>(90% max)</td><td>Maximum value to which image is
% scaled</td></tr>
% <tr><td>PD</td><td>0</td><td>Declares that input image is positive
% definite.</td></tr>
% <tr><td>cbmode</td><td>0</td><td>Specifies whether to use custom colorbar
% axis labels.</td></tr>
% <tr><td>clabels</td><td>([-90% max,  90% max])</td><td>Colorbar axis
% labels. When cbmode==1, min defaults to 0 if PD==1, both default to +/-
% Scale if supplied. When cbmode==0, then cblabels dictates colorbar axis
% limits.</td></tr>
% <tr><td>cbticks</td><td>(none)</td><td>When cbmode==1, specifies
% positions of tick marks on colorbar axis.</td></tr>
% <tr><td>slices</td><td>(center frames)</td><td>Select which slices are
% displayed. If empty, activates interactive navigation.</td></tr>
% <tr><td>slices_type</td><td>'idx'</td><td>Use MATLAB indexing ('idx') for
% slices, or spatial coordinates ('coord') as provided by
% invoVol.</td></tr>
% <tr><td>orientation</td><td>'t'</td><td>Select orientation of volume.
% 't' for transverse, 's' for sagittal.</td></tr>
% <tr><td>kernel</td><td>[1]</td><td>A sampling kernel for the time trace
% plot. Other options: 'gaussian' | 'cube' | 'sphere'.</td></tr>
% </table></html>
%
% Note: |applycmap| has further options for using |params| to specify
% parameters for the fusion, scaling, and colormapping process.
%
%% Dependencies
% <applycmap_help.html applycmap>
% 
%% See Also
% <PlotInterpSurfMesh_help.html PlotInterpSurfMesh> | <PlotSlices_help.html
% PlotSlices>

