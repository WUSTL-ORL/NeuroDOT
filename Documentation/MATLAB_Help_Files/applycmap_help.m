%% applycmap
% Performs color mapping and fuses images with anatomical models.
%
%% Description
% |img_out = applycmap(img_in)| fuses the N-D image array |img_in| with a
% default flat gray background and applies a number of other default
% settings to create a scaled and colormapped N-D x 3 RGB image array
% |img_out|.
%
% |img_out = applycmap(img_in, anat)| fuses the image array |img_in| with
% the anatomical atlas volume input |anat| as the background.
%
% |img_out = applycmap(img_in, anat, params)| allows the user to specify
% parameters for plot creation.
%
%% Visualization Parameters
% "params" fields that apply to this function (and their defaults):
% 
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>TC</td><td>0</td><td>Direct map integer data values to defined
% color map ("True Color").</td></tr>
% <tr><td>DR</td><td>1000</td><td>Dynamic range.</td></tr>
% <tr><td>Scale</td><td>0.9</td><td>Maximum value, as a decimal percentage,
% to which image is scaled.</td></tr>
% <tr><td>PD</td><td>0</td><td>Declares that input image is positive
% definite.</td></tr>
% <tr><td>Cmap</td><td></td><td>Color maps.</td></tr>
% <tr><td>.P</td><td>'jet'</td><td>Colormap for positive data
% values.</td></tr>
% <tr><td>.N</td><td>(none)</td><td>Colormap for negative data
% values.</td></tr>
% <tr><td>.flipP</td><td>0</td><td>Logical, flips the positive
% colormap.</td></tr>
% <tr><td>.flipN</td><td>0</td><td>Logical, flips the negative
% colormap.</td></tr>
% <tr><td>Th</td><td></td><td>Thresholds.</td></tr>
% <tr><td>.P</td><td>(25% max)</td><td>Value of min threshold to display
% positive data values.</td></tr>
% <tr><td>.N</td><td>(-Th.P)</td><td>Value of max threshold to display
% negative data values.</td></tr>
% <tr><td>BG</td><td>[0.5, 0.5, 0.5]</td><td>Background color, as an RGB
% triplet.</td></tr>
% <tr><td>Saturation</td><td>(none)</td><td>Field the size of data with
% values to set the coloring saturation. Must be within range [0,
% 1].</td></tr>
% </table></html>
% 
%% See Also
% <PlotSlices_help.html PlotSlices> | <PlotInterpSurfMesh_help.html
% PlotInterpSurfMesh> | <PlotLRMeshes_help.html PlotLRMeshes> |
% <PlotCapMeanLL_help.html PlotCapMeanLL>