function varargout = Cap_Fitter(varargin)

% CAP_FITTER An interactive tool for fitting optode grids onto the skull.
% 
%   CAP_FITTER is an interactive GUI that allows the user to load skull
%   mesh and cap grid files, re-orient them with an affine transformation,
%   and then run a spring-fit algorithm on them. The final result will then
%   be saved in a file containing the original "info" structure, the
%   resulting one, the mesh they were fit onto, and an "ops" structure
%   recording the operations performed.
% 
% Dependencies: PLOTMESHSURFACE, PLOTCAP, SPRINGFIT_OPTODES,
% OPERATE_OPTODES, SCALE_CAP, ROTATE_CAP, ROTATION_MATRIX.% 
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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Cap_Fitter_OpeningFcn, ...
    'gui_OutputFcn',  @Cap_Fitter_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Cap_Fitter is made visible.
function Cap_Fitter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cap_Fitter (see VARARGIN)

% Choose default command line output for Cap_Fitter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Cap_Fitter wait for user response (see UIRESUME)
% uiwait(handles.Cap_Fitter);

% Add degree symbols.
handles.RotXText.String = ['x', char(176)];
handles.RotYText.String = ['y', char(176)];
handles.RotZText.String = ['z', char(176)];

global info info2 mesh ops pathname

pathname = pwd;
info = [];
info2
mesh = [];
ops = [];

clc


% --- Outputs from this function are returned to the command line.
function varargout = Cap_Fitter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadHeadMesh_Callback(hObject, eventdata, handles)

global mesh pathname

try
    % Open file dialog.
    [filename, pathname, idx] = uigetfile('*.mat', 'Select Head Mesh File', pathname);
    assert(logical(idx))
    
    % Load file.
    S = load(fullfile(pathname, filename));
    if isfield(S, 'mesh')
        mesh = S.mesh;
    else
        error('*** Error: No mesh data detected. Please try another file. ***')
    end
    
    % Display the mesh in flat mode.
    params.Nregs = 0;
    params.fig_handle = gcf;
    params.BG = 'w';
    PlotMeshSurface(mesh, [], params)
    view(163, -86)
    handles.Az.String = num2str(163);
    handles.El.String = num2str(-86);
    AssignRotateCallback
    
    % Update the rest of the panel.
    EnableNav(handles)
    handles.LoadCapGrid.Enable = 'on';
    handles.HeadText.String = filename;
catch err
    pathname = [];
%     rethrow(err)
end


function LoadCapGrid_Callback(hObject, eventdata, handles)

global info pathname

try
    % Open file dialog.
    [filename, pathname, idx] = uigetfile('*.mat', 'Select Cap Grid File', pathname);
    assert(logical(idx))
    
    % Load file.
    S = load(fullfile(pathname, filename));
    
    % Detect data type.
    if isfield(S, 'info')
        info = S.info;
    elseif isfield(S,'grid')
        info.optodes.spos3=S.grid.spos3;
        info.optodes.dpos3=S.grid.dpos3;
        info.pairs.Src=1:S.grid.srcnum;
        info.pairs.Det=1:S.grid.detnum;
        info.optodes.plot3orientation.i='R2L';
        info.optodes.plot3orientation.j='P2A';
        info.optodes.plot3orientation.k='D2V';
    else
        
%         error('*** Error: No grid data detected. Please try another file. ***')
    end
    
    % Display the grid.
    ResetTransform_Callback(hObject, eventdata, handles)
    
    % Update the rest of the panel.
    EnablePanel(handles);
    handles.CapText.String = filename;
catch err
    pathname = [];
%     rethrow(err)
end


function UpdateCap_Callback(hObject, eventdata, handles)

UpdateData(handles)

UpdateAxes

handles.RunSpringFit.Enable = 'on';


function ResetCapView_Callback(hObject, eventdata, handles)

view(163, -86)


function ResetTransform_Callback(hObject, eventdata, handles)

handles.RotX.String = '0';
handles.RotY.String = '0';
handles.RotZ.String = '0';

handles.TransX.String = '0';
handles.TransY.String = '0';
handles.TransZ.String = '0';

handles.FlipX.Value = 0;
handles.FlipY.Value = 0;
handles.FlipZ.Value = 0;

handles.ScaleM.String = '1';

UpdateData(handles)
UpdateAxes


function RunSpringFit_Callback(hObject, eventdata, handles)

global info2 mesh

Ns = size(info2.optodes.spos3, 1);
Nd = size(info2.optodes.dpos3, 1);

% Convert current transformed data to ND2 format (COMPATIBILITY FOR
% SPRINGFIT).
info2 = converter_info(info2, 'ND2 to ND1');
info2.iRad.r = info2.iRad.ir;
info2.iRad.meas = [];
info2.iRad.meas(:, 1) = repmat([1:Ns]', Nd, 1);
info2.iRad.meas(:, 2) = floor(1:1/Ns:Nd+(Ns-1)/Ns)';

% Load an old ND1 rad file.
rad = load('radius_Adult_96x92.mat');

% Run the spring fit on transformed data.
info2.tpos = springfit_optodes(mesh, rad, info2.grid.spos3,...
    info2.grid.dpos3, [], []);

% % This is only for development purposes.
% load('springfit_SAVE_STATE.mat')
% close(gcf)

info2 = converter_info(info2, 'ND1 to ND2');
info2.optodes.spos3 = info2.misc.tpos(1:96, :);
info2.optodes.dpos3 = info2.misc.tpos(97:end, :);

% Update the axes.
UpdateAxes


function SaveCapPosition_Callback(hObject, eventdata, handles)

global info info2 mesh ops pathname %#ok<NUSED>

try
    % Save file dialog.
    [filename, pathname, idx] = uiputfile('*.mat', 'Save Cap Fit Results', pathname);
    assert(logical(idx))
    
    infos.info_in = info;
    infos.info_out = info2;
    infos.info_out.optodes = rmfield(infos.info_out.optodes, 'tpos');
    
    % Save file.
    save(fullfile(pathname, filename), 'infos', 'mesh', 'ops')
    
catch err
    rethrow(err)
end


function AzL1_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az + 1, el)
handles.Az.String = num2str(az + 1);


function AzR1_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az - 1, el)
handles.Az.String = num2str(az - 1);


function ElL1_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az, el + 1)
handles.El.String = num2str(el + 1);


function ElR1_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az, el - 1)
handles.El.String = num2str(el - 1);


function AzL10_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az + 10, el)
handles.Az.String = num2str(az + 10);


function AzR10_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az - 10, el)
handles.Az.String = num2str(az - 10);


function Az_Callback(hObject, eventdata, handles)

[az, el] = view;
view(str2double(hObject.String), el)


function ElR10_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az, el + 10)
handles.El.String = num2str(el + 10);


function ElL10_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az, el - 10)
handles.El.String = num2str(el - 10);


function El_Callback(hObject, eventdata, handles)

[az, el] = view;
view(az, str2double(hObject.String))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UpdateAxes

global mesh info2

% Update the plot.
cla

params.Nregs = 0;
params.fig_handle = gcf;
PlotMeshSurface(mesh, [], params)

hold on % KEEP THIS HERE!

params = [];
params.fig_handle = gcf;
params.dimension = '3D';
params.BG = 'w';
PlotCap(info2, params)
AssignRotateCallback


function AssignRotateCallback

h = rotate3d(gcf);
h.ActionPostCallback = @ActionPost;


function UpdateData(handles)

global ops info info2

% Collect values.
ops.rot_x = str2double(handles.RotX.String);
ops.rot_y = str2double(handles.RotY.String);
ops.rot_z = str2double(handles.RotZ.String);

ops.trans_x = str2double(handles.TransX.String);
ops.trans_y = str2double(handles.TransY.String);
ops.trans_z = str2double(handles.TransZ.String);

ops.flip_x = handles.FlipX.Value;
ops.flip_y = handles.FlipY.Value;
ops.flip_z = handles.FlipZ.Value;

ops.scale_m = str2double(handles.ScaleM.String);

info2 = operate_optodes(info, ops);


function EnableNav(handles)

handles.AzText.Enable = 'on';
handles.Az.Enable = 'on';
handles.AzL1.Enable = 'on';
handles.AzR1.Enable = 'on';
handles.AzL10.Enable = 'on';
handles.AzR10.Enable = 'on';
handles.ElText.Enable = 'on';
handles.El.Enable = 'on';
handles.ElL1.Enable = 'on';
handles.ElR1.Enable = 'on';
handles.ElL10.Enable = 'on';
handles.ElR10.Enable = 'on';


function EnablePanel(handles)

handles.RotText.Enable = 'on';
handles.RotX.Enable = 'on';
handles.RotY.Enable = 'on';
handles.RotZ.Enable = 'on';
handles.RotXText.Enable = 'on';
handles.RotYText.Enable = 'on';
handles.RotZText.Enable = 'on';

handles.TransText.Enable = 'on';
handles.TransX.Enable = 'on';
handles.TransY.Enable = 'on';
handles.TransZ.Enable = 'on';
handles.TransXText.Enable = 'on';
handles.TransYText.Enable = 'on';
handles.TransZText.Enable = 'on';

handles.FlipText.Enable = 'on';
handles.FlipX.Enable = 'on';
handles.FlipY.Enable = 'on';
handles.FlipZ.Enable = 'on';

handles.ScaleText.Enable = 'on';
handles.ScaleM.Enable = 'on';
handles.ScaleMText.Enable = 'on';

handles.UpdateCap.Enable = 'on';
handles.ResetCapView.Enable = 'on';
handles.ResetTransform.Enable = 'on';
handles.SaveCapPosition.Enable = 'on';


function ActionPost(hObject, event)

handles = guihandles(hObject);

[az, el] = view;
handles.Az.String = num2str(az);
handles.El.String = num2str(el);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unused Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FlipX_Callback(hObject, eventdata, handles)
function FlipY_Callback(hObject, eventdata, handles)
function FlipZ_Callback(hObject, eventdata, handles)
function RotX_Callback(hObject, eventdata, handles)
function RotY_Callback(hObject, eventdata, handles)
function RotZ_Callback(hObject, eventdata, handles)
function TransX_Callback(hObject, eventdata, handles)
function TransY_Callback(hObject, eventdata, handles)
function TransZ_Callback(hObject, eventdata, handles)
function ScaleM_Callback(hObject, eventdata, handles)
function RotX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RotY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RotZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TransX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TransY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TransZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ScaleM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Az_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function El_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
