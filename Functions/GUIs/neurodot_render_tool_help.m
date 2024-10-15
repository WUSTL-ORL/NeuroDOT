%NEURODOT_RENDER_TOOL
%
% This utility allows you to interactively select a nifti or 4dfp overlay 
% and render the data on a default or custom cortical mesh. Alternatively,
% data may be passed to the app directly from the command line or from
% within a script. See RENDER TOOL STARTUP (below) for parameter options.
%
% Note NeuroDOT must be installed and included in your Matlab path to use
% the render tool.
%
% The following GUI controls are provided:
%
%   Overlay / Select new -- select the overlay to display. Recognized file
%   formats are nifti (.nii), 4dpf (.img/.ifh), and native NeuroDOT (.mat).
%   If a 4D file is selected, a frame selection control is added to the GUI
%
%   Mesh / Select new -- select an alternate cortical mesh to display,
%   such as a mesh created from a subjects native structural image. 
%   The selected file must be a .mat file containing valid NeuroDOT
%   cortical mesh data in structures named MESH_L and MESH_R, or MNIl
%   and MNIr, or Anat.CtxL and Anat.CtxR (all are case sensitive)
%
%   NB: The mesh and overlay must be in the same space. A mismatch may
%   cause the overlay to be rendered incorrectly. The least error-prone
%   approach is to use overlays normalized to an MNI template, as that 
%   is usually compatible with the default cortical mesh
%
%   Color Scale -- sets the scale range of the colormap. Smaller values 
%   bring out hotspots (increases saturation of active sites); larger 
%   values decrease contrast. Absolute or relative scaling are selected 
%   using the radio buttons (in relative scaling, the specified scale 
%   value is multiplied by the maximum value of the current overlay)
%
%   You may change scaling using either the scale slider or by 
%   entering a new value in the textbox to the right of the slider. If 
%   the new value exeedes the current maximum value of the slider, it
%   becomes the new maximum scale value
%
%   (-) Thresh -- set the negative threshold, as a percent of the maximum
%   negative value of the color scale. As the threshold is increased, 
%   less data in the overlay is displayed
%
%   (+) Thresh -- set the positive threshold
%
%   NB: Thresholding is relative to the color scale. Activation may still 
%   be apparent even at 100% threshold if the overlay contains values that 
%   exceed the maximum and/or minimum of the current color scale range
%
%   Colorbar [on/off] -- determines whether a colorbar is displayed
%
%   Colormap [various] -- sets the colormap
%
%   View [lat(eral)/dorsal/post(erior)] -- sets the view orientation
%
%   Inflation [none/some/max] -- sets the cortex template inflation
%
%   NB: If the active cortical mesh does not contain inflated meshes, the 
%   Inflation menu will be unavailable
%
%   Overlap interpolation [linear/nearest] -- sets overlap interpolation
%   used to display the overlay
%
%   Only plot positive overlay values -- selecting this option will ignore
%   any negative values in the overlay. The colormap will be rescaled to
%   [ 0 max(overlay>0) ]. In NeuroDOT documentation, this is refered to as
%   the overlay being "positive definitive"
%
%   Print current settings to command window -- prints the "params" struct,
%   the positive and negative thresholds, and overlay and template info
%   to the command window to document settings. This can be convenient
%   to identify parameter values useful in subsequent scripting
%
%   Restore Defaults -- restore all setting to internal default values
%
%   Help -- display this help text in a scrollable dialog window. Note the
%   help window can be moved and resized as needed
%
%   Quit -- quit the app (closes the render and help windows, if open)
%
% NB: Menus in the render window may be used to save the image, zoom, etc.
%
%
% RENDER TOOL STARTUP
%
% In addition to double-clicking the app icon, the rendertool may be
% started from the command line or from within an mfile. The following 
% calling conventions are implemented:
%
%   1) neurodot_render_tool (no parameters). The tool opens as if the icon
%   was double-clicked. The overlay, mesh, etc are then selected using the 
%   GUI
%
%   2) neurodot_render_tool(overlay_filename). The tool renders the
%   specified overlay using the default settings (which can then be changed
%   interactively). The file must exist or the app exits with an error
%   message
%
%   3) neurodot_render_tool(overlay_data,overlay_header). Similar to option
%   #2, except the overlay data and header struct are passed instead of a
%   filename. Rudimentary input checking is applied and the app will exit
%   with an error message if improper data is detected
%
%   4) neurodot_render_tool(overlay_data,overlay_header,lmesh,rmesh).
%   Similar to option #3, except structures for the left and right
%   cortical meshes are also passed
%
% Additionally, the final (or only) parameter passed can be a structure
% defining startup values for rendering settings (scaling, colormap, etc).
% This defines four additional possible startup commands:
%
%  >> rendertool(params);
%  >> rendertool(overlay_filename,params);
%  >> rendertool(overlay_data,overlay_header,params);
%  >> neurodot_render_tool(overlay_data,overlay_header,lmesh,rmesh,params)
%
% Two versions of the parameter struct are accepted. The first defines a
% set of (hopefully) self-documenting fieldnames:
%
%     colormap                      = [ any matlab colormap except "vga" ]
%     color_scale                   = [ float ]
%     color_scale_mode              = [ 'absolute' | 'relative' ]
%     positive_threshold_percent    = [ integer: 0-100 ]
%     negative_threshold_percent    = [ integer: 0-100]
%     show_colorbar                 = [ true | false ]
%     view                          = [ 'lat' | 'dorsal' | 'post' ]
%     cortex_inflation              = ['none' | 'some? | 'max' ]
%     interpolation                 = [ 'linear' | 'nearest' ]
%     plot_positive_only            = [ true | false ]
%
% Example usage:
%
%     >> params = [ ];
%     >> params.colormap = 'parula';
%     >> params.color_scale = 0.5;
%     >> neurodot_render_tool('NeuroDOT_Render_Tool_Example.nii',params);
%
% Alternatively, the NeuroDOT internal rendering parameter structure can
% be passed. For example: 
%
%     >> params = [ ];
%     >> params.Scale: 0.0375
%     >> params.Cmap = 'jet';
%     >> params.view = 'lat'
%     >> params.CBar_on = 1;
%     >. params.PD = 0;
%     >> params.ctx = 'std';
%     >> params.OL = 0;
%     >> params.Th.P: 0.01;
%     >> params.Th.N: -0.01;
%     >> neurodot_render_tool('NeuroDOT_Render_Tool_Example.nii',params);
%
% See PlotInterpSurfMesh.m for a description of these fields.
%
% PARAMETER STRUCTURE NOTES
%
% - If the passed structure does not define a parameter, the default 
% setting will be used
%
% - Although some sanity checking is performed on passed parameters,
% it's possible to specify invalid settings that crash the render tool
%
% - If you specify a colormap that is not included in the default set, 
% it will be added to the colormap pulldown menu
%
%
% KNOWN ISSUES
%
% - positive and negative thresholding interact with each other in
% some colormaps (changing one alters the other)
%
% - as of this writing, file type filtering is broken in the Matlab file
% selection dialog box under OS-X. As such, no filtering is used, which 
% potentially allows selection of inappropriate file types (e.g., a 
% non-mat file when selecting a custom cortical mesh). The file type is
% instead checked after selection and an error message is displayed if an
% incorrect type is detected
%
% - selecting a new overlay or mesh file when the help window is open may
% cause two file selection dialog boxes to appear. Use one to select the 
% file, then select "Cancel" to close the other
%
% - NeuroDOT may incorrectly warn "The Overlay has only elements equal
% to zero" when rendering some overlays
%

% CHANGE HISTORY
%
% 2024 [MSJ] - new
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
