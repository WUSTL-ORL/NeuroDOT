%% Script for Participant Specific Head Modeling
% This script follows the tutorial:
% "Tutorial_for_Participant_Specific_Head_Modeling.pptx"

% This script utilizes FreeSurfer 7.2:
% https://surfer.nmr.mgh.harvard.edu/fswiki/rel7downloads
% DOI: 10.1016/j.neuroimage.2012.01.021 (Fischl 2012)

% Inputs (Section 0):
%   dataroot: Path to the folder containing the participant's MRI data
%   fn: Filename of the T1 MRI 
%   t2Mode: Flag for using a T2 MRI in addition to the T1 ('on' or 'off')
%   fn_T2: Filename of the T2 MRI (optional)
%       NOTE: The T2 MRI must be in the same orientation and in exact
%       alignment with the T1 MRI.
%   padname: the name of the imaging array used which corresponds to the
%   correct "pad file"
%   mode: 'participant' or 'atlas' options for head modeling procedure
%   pt: participant ID
%   atlasDir: Path to FreeSurfer output directory
%   codeDir: Path to location of fsLR code

% This script performs the following tasks:
%   1. Head segmentation and data registration with FreeSurfer
%       1a. Run FreeSurfer and fsLR and unpack their outputs
%       1b. Segment the head with NeuroDOT
%   2. Computational head mesh generation and array alignment
%       2a. Make a low-density head mesh
%       2b. Use AlignMe to align the array Make a NIRFAST-compliant with the mesh
%       2c. Relax optodes on the head mesh and save relaxed optode positions
%       2d. Optimize mask for a targeted computational mesh
%       2e. Make a NIRFAST-compliant high-density head mesh
%       2f. Relax optodes on the NIRFAST-compliant mesh for light modeling
%   3. Light Modeling
%       3a. Make “A” sensitivity matrix with NIRFASTer
%       3b. Package metadata and save with “A”
%       3c. Update bookkeeping in “A” to save space and align with system data
%       3d. Visualize sensitivity profile and flat field reconstruction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0: Parameters and initialization
dataroot = '';           % Make a new directory and set your path here

dataout = ''             % Select the directory for fsLR outputs to be placed. 
                         % Note: dataout must be a full path,
                         % i.e. C:/Users/username/dataroot/dataout.
                         % Do NOT place a "/" at the end of the path 
                                                            
fn = '';                 % The T1-weighted MRI sent to FreeSurfer  (omit file extension):
file_type = 'nii';       % File extension of T1 (nii or 4dfp)
T1_path = '';            % Hardcoded, do not change

% If you have a T2-weighted MRI which is in register with the T1 MRI, 
% the name (with no extensions): 

t2Mode = 'off';
if strcmp(t2Mode, 'on')
    fn_T2 = ''; 
end

padname = 'LUMO_adult_visual';         %  the imaging array to use
mode = 'participant';                  % 'atlas' for atlas-based head model,'participant' for participant head model
pt = '';                               % Participant ID here (NOTE: this must be changed each time the script is run)

% Path to fsLR script
fsLRDir = 'NeuroDOT/Functions/Light_Modeling/Participant_Specfic_Head_Modeling_fs_LR_Script.sh'; % Hardcoded, do not change

% Path to Workbench command 
workbenchDirectory = '/usr/local/workbench/bin_rh_linux64'; % Set the path to your installation of workbench here

% Volumetric visualization parameters
p.Cmap = 'jet'; p.Scale = 1; p.Th.P = 0.5; p.PD = 0;

% Mesh visualization parameters
pM = p; pM.Scale = 5; pM.Cmap = 'gray'; pM.BG = [0,0,0]; pM.PD = 1; pM.EdgeColor = 'none';

% Sensitivity profile visualization parameters
pA.PD = 1; pA.Scale = 2; pA.Th.P = 1e-2; pA.Th.N = -pA.Th.P;


%% Section 1
%% 1a: Volumetric Segmentation with FreeSurfer and realignment of data

% FREESURFER SETUP -----------------------------------------------------------

% full path to top level install directory
FSsetup.freesurferdir = '/usr/local/freesurfer';

% shell used to run FreeSurfer
FSsetup.freesurfershell = 'bash';

% full path to FreeSurfer setup script, run before any FS command
FSsetup.freesurfersetup = '/usr/local/freesurfer/SetUpFreeSurfer.sh';

% full path to FreeSurfer environmental setup script, run before any FS command
FSsetup.freesurferenvironment = '/usr/local/freesurfer/FreeSurferEnv.sh';

% full path to FreeSurfer "subjects" directory - CHECK: now unused?
% maybe a user will put these files somewhere else??
FSsetup.freesurfersubjectsdir = '/usr/local/freesurfer/subjects';

% freesurfer output format ('nii' or 'nii.gz')
FSsetup.freesurferoutputformat = 'nii.gz';

pathtostructural = fullfile(dataroot,[fn '.nii']);
outputdirectoryname = pt; % created by freesurfer, should not exist
parentoutputdirectory = dataroot;


FScmd = sprintf('recon-all -all -i %s -s %s -sd %s', pathtostructural, outputdirectoryname, parentoutputdirectory);

t=tic;          
run_freesurfer_command(FScmd,FSsetup);
toc(t)


%% Unpack FreeSurfer Outputs
% Put T1, brainmask, and aseg on mpr1
cd([dataroot,'/freesurfer/','/mri/'])
[status,response]=system(['mri_vol2vol --mov T1.mgz --targ ',...
    'rawavg.mgz --regheader --o T1-in-rawavg.mgz']);
[status,response]=system(['mri_label2vol --seg aseg.mgz --temp ',...
    'rawavg.mgz --o aseg-in-rawavg.mgz --regheader aseg.mgz']);  %#ok<*ASGLU>
[status,response]=system(['mri_label2vol --seg brainmask.mgz --temp ',...
    'rawavg.mgz --o brainmask-in-rawavg.mgz --regheader brainmask.mgz']);

% MGZ 2 NIFTI
disp('Convert *.mgz files to *.nii');          % Convert mgz to nifti
[status,response]=system(['mri_convert T1-in-rawavg.mgz T1.nii']);
[status,response]=system(['mri_convert aseg-in-rawavg.mgz aseg.nii']);
[status,response]=system(['mri_convert brainmask-in-rawavg.mgz brainmask.nii']);

% Now load in volumes to complete segmentation: T1, aseg, brainmask, t2w
[T1,infoT1] = LoadVolumetricData('T1', [],'nii'); % Output from FreeSurfer 
T1=T1./max(T1(:));
[bm,info_bm] = LoadVolumetricData(['brainmask'], [],'nii'); % Load Brain mask
[aseg, info_aseg] = LoadVolumetricData(['aseg'], [],'nii'); % Load Anatomical segmentation volume

% Navigate back to original T1 data directory
cd(dataroot)
[mpr1,infompr1] = LoadVolumetricData(fn, [],'nii'); % Load original T1 input file 

% This section will be run if you are using a T2 in addition to a T1
switch t2Mode
    case 'on'
        % Load T2 If T2 file is present and aligned with T1
        [T2_on_T1,infoT2] = LoadVolumetricData(fn_T2, [],'nii'); 
        T2_on_T1 = T2_on_T1./max(T2_on_T1(:));
        % If the T1 and T2 are in register, the CSF of the T2 will fit exquisitely
        % in the correct spots of the T1.
        PlotSlices(T1, infoT1, T2_on_T1, p) 
end

% Visualize each to ensure they are co-registered
PlotSlices(T1,infoT1)
PlotSlices(aseg,infoT1)
PlotSlices(T1,infoT1,p,bm)


%% Run fsLR to put left and right hemispheres in correspondence
cd(dataroot)
t1=tic;
status= evalc(['system([''/',fsLRDir,' ',dataroot,' ',pt,' ',dataroot,'/',pt,...
    '/mri ',dataroot,'/standard_mesh_atlases ', dataroot, '/freesurfer/mri ', ...
    dataroot,'/freesurfer',' ',dataroot,'/',pt,' ',dataroot,'/',pt,'/mri/T1.nii ',...
    dataroot,'/standard_mesh_atlases ', '164', ' ','32', ' ',...
    workbenchDirectory,'''])']);
disp('<<< Pipeline to  fs_LR')
toc(t1)


%%
evalc(['system([''mkdir ',dataout,'/fs_LR/''])']);
resArray = [32, 164];
for i = 1:2
    res = resArray(i);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.L.pial.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.L.inflated.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/', pt,'.L.very_inflated.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.L.flat.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.L.curvature.', num2str(res),'k_fs_LR.shape.gii',' /',...
        dataout,'/fs_LR/''])']);

    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.R.pial.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.R.inflated.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/', pt,'.R.very_inflated.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.R.flat.', num2str(res),'k_fs_LR.surf.gii',' /',...
        dataout,'/fs_LR/''])']);
    sys = evalc(['system([''cp ',dataroot,'/standard_mesh_atlases/','fsaverage_LR', num2str(res),'k/',pt,'.R.curvature.', num2str(res),'k_fs_LR.shape.gii',' /',...
        dataout,'/fs_LR/''])']);
end

subjSessID = ['sub-',pt];
subjsSessID = pt;

% Create Anat structure
params.res = 'high'; %set to 'low' for 32k or 'high' for 164k
switch params.res
    case 'low'
        Anat=Convert_FSLR_WB_Ctx2mat(pt, [dataout, '/fs_LR'],params);
    case 'high'
        Anat=Convert_FSLR_WB_Ctx2mat(pt, [dataout, '/fs_LR'],params);
end
cd([dataout,'/fs_LR'])

% Visualize pial, inflated, very inflated, and flat surfaces
res = '164k'; % 164k or 32k 
load([pt,'_', res, '_ctx.mat']);

% Initialize fsLR visualization parameter structure
params_fsLR = struct;

% View pial surface
params_fsLR.ctx = 'std';
PlotLRMeshes(Anat.CtxL,Anat.CtxR,params_fsLR)

% View inflated surface
params_fsLR.ctx = 'inf';
PlotLRMeshes(Anat.CtxL,Anat.CtxR,params_fsLR)

% View very inflated surface
params_fsLR.ctx = 'vinf';
PlotLRMeshes(Anat.CtxL,Anat.CtxR,params_fsLR)

% View flat surface
params_fsLR.ctx = 'flat'; params_fsLR.reg = 1;
params_fsLR.reg = 1; params_fsLR.FaceColor = 'interp'; params_fsLR.EdgeColor = 'none';
params_fsLR.view = 'lat';
PlotMeshSurface(Anat.CtxL,params_fsLR); 
set(gca, 'color','k', 'Xcolor','k', 'YColor','k'); set(gcf, 'color','k');
PlotMeshSurface(Anat.CtxR,params_fsLR)
set(gca, 'color','k', 'Xcolor','k', 'YColor','k'); set(gcf, 'color','k');


%% 1b: Segment the rest of the head and make the whole head mask
tag=[''];
params.killY=[1:10]; % Voxels to kill on neck to save space
params.head_bg_seg_type='T1';
params.Skull_th=0.15;
params.Head_th = 0.15;
switch t2Mode
    case 'on'
        [mask,params]=Segment5R_fs(T1,T2_on_T1,bm,aseg,infoT1,params); % Only if you have T2
    case 'off'
        [mask] = Segment5R_fs_noT2(T1, bm, aseg, infoT1, params);
end

pM1 = pM;
pM1.Cmap = 'jet';
PlotSlices(mask,infoT1,pM1)

% Save Segmented Volume as NIFTI file
cd(dataroot)
SaveVolumetricData(mask, infoT1,[pt,'_Segmented_Mask_',tag],[],'nii');

% Clear up workspace except for variables needed for Light Modeling
clearvars -except dataroot dataout tag padname pt mode res p pM pM1 pA
close all


%% Section 2: Load a Segmented Volume and generate head mesh
[mask, infoT1] = LoadVolumetricData([pt,'_Segmented_Mask_',tag],[],'nii');
PlotSlices(mask,infoT1,pM1)       % Visualize the segmented mask: 
        % Note, PlotSlices is an interactive plot. 
        % Hit q or the middle mouse button to quit


%% 2a. Generate low density head mesh
% Parameters for generating your mesh
meshname=['LD_Mesh'];      % Provide a name for your mesh name here
param.facet_distance=5;    % Node position error tolerance at boundary
param.facet_size=3;     % boundary element size parameter
param.cell_size=5;      % Volume element size parameter
param.info=infoT1;
param.Mode=1;                    % make simple mesh with no region labels

switch mode
    case 'atlas'
        param.Offset = [0,0,0]; 
        tic;meshLD=NirfastMesh_Region(mask,meshname,param);toc
        % Put coordinates back in true space (Only run if mesh coordinates
        % are not centered at zero)
        meshLD.nodes=change_space_coords(meshLD.nodes,infoT1,'coord'); 
    case 'participant'
        param.Offset=infoT1.center; 
        tic;meshLD=NirfastMesh_Region(mask,meshname,param);toc
end
PlotMeshSurface(meshLD,pM) %visualize in index (non-negative) space
view([135,30]) %view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is

% Look at the inside
m1.nodes=meshLD.nodes;m1.elements=meshLD.elements;
m2=CutMesh(m1,find(m1.nodes(:,2)>70)); % Cut LD mesh in half
[Ia,Ib]=ismember(m2.nodes,meshLD.nodes,'rows');Ib(Ib==0)=[];
PlotMeshSurface(m2,pM);               % Visualize back half of LD mesh
% View the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is
%now the mesh should be centered at (0,0,0).
% view([135,30])                        
view([0,0]); % Posterior view


%% 2b. Load initial Pad file depending on data structure:
load(['Pad_',padname,'.mat']); % Pad file

tpos=cat(1,info.optodes.spos3,info.optodes.dpos3);
rad=info;
% % Array visualization parameters
Ns=size(info.optodes.spos3,1);
Nd=size(info.optodes.dpos3,1);
paramsFoci.color = cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1));
paramsFoci.color(1,:) = [1 0.4 0.6]; % Pink for S1
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; % Light blue for D1
paramsFoci.radius = 2; % Set radius of spheres to 2

% Visualize Pad
figure;Draw_Foci_191203(tpos, paramsFoci);view([-40,30]) %3D plot of pad
params_cap.eeg_style = 1; %turn on EEG style for this cap
PlotCap(info, params_cap) %2D plot of pad

% Visualize LD mesh with optode positions 
% Optodes have not been relaxed onto mesh yet
% % Array visualization parameters
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tpos, paramsFoci);


%% Make sure only NN1 and NN2 exist for meas < 6cm
% Visualize scatterplot of Nearest neighbors for all measurements under 6cm
% Set cutoff at 6cm (60mm), we won't use any measurements larger than this
keep = info.pairs.r3d < 60; 
figure;
scatter(info.pairs.NN(keep),info.pairs.r3d(keep));
xlabel('NN');ylabel('r3d');title('NN vs SD separation');

% If the plot shows measurements for NN3 or greater, run the following 
% lines (See slide 12 of tutorial powerpoint)
% Set all meas with separation < 3cm to NN1 and all meas in between 3cm
% and 6cm separation to NN2
info.pairs.NN(info.pairs.r3d < 30) = 1; 
info.pairs.NN(info.pairs.r3d >= 30 & info.pairs.r3d < 60) = 2; 
figure; % Visualize
scatter(info.pairs.NN(keep),info.pairs.r3d(keep));
xlabel('NN');ylabel('r3d');title('NN vs SD separation');


%% 2b-c. AlignMe Section:  Move grid from arbitrary location to approximate target on mesh, Relax grid on head and view
atlasFiducials = [- 0.65,  -84.1, -31.88; ... % Nasion --> EEGPts(1,:)
                  - 0.65, 117.69, -11.78; ... % Inion  EEGPts(329,:)
                   80.78,  15.95, -41.89; ... % LPA    EEGPts(155,:)
                  -80.78,  15.95, -41.89; ... % RPA    EEGPts(175,:)
                   0.233,   9.63, 97.296];    % Cz     EEGPts(165,:)
               
% Create an instance of our custom DataStorage HANDLE class to store variables
ds = DataStorage(); 

% Input structure
ds.dI.tpos = tpos;       
ds.dI.mesh = meshLD;
ds.dI.pad = info;
ds.dI.pM = pM;
ds.dI.paramsFoci = paramsFoci;
ds.dI.Ns = Ns;
ds.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application, passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
close all
myapp=AlignMe_2020b(ds);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew = ds.dO.tpos2_relaxed;

% visualize mesh with relaxed optodes
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew,paramsFoci)
view(0,0) %posterior view


%% 2d. Generate smaller mask based on optode locations
maskCrop=IceCreamScoopCap4Mask(infoT1,tposNew,mask,40); %40 = mm of mesh to scoop out
PlotSlices(maskCrop,infoT1,pM1) %visualize smaller mask


%% 2e. Generate High Density Head Mesh (Can take up to 15 min)
%If you get an error when running NirfastMesh_Region
%try changing the mesh name and clearing your output directory of all files
% Provide a name for your mesh name here, if making multiple meshes,
% provide a different name for each mesh
meshname=[''];      
% Node position error tolerance at boundary
param.facet_distance=2.0;   
% boundary element size parameter
param.facet_size=0.8;       
% Volume element size parameter - 
%(for this mesh, 1.2 is the lowest you can set the cell_size for to 
% have less than 999,999 nodes)
param.cell_size=1.5;        
param.Mode=0;               % Set Mode=0 to make nirfast compliant mesh
param.info=infoT1;          % Make sure info is infoT1 like with LD mesh
switch mode
    case 'atlas'
        param.Offset = [0,0,0]; 
        tic;meshHD=NirfastMesh_Region(maskCrop,meshname,param);toc
        % Put coordinates back in true space
        meshHD.nodes=change_space_coords(meshHD.nodes,infoT1,'coord'); %for actual mesh
        visMeshHD.nodes = meshHD.nodes; %make copy of mesh that only contains nodes and elements for visualization purposes
        visMeshHD.elements = meshHD.elements;
    case 'participant'
        param.Offset=infoT1.center; 
        tic;meshHD=NirfastMesh_Region(maskCrop,meshname,param);toc     %for actual mesh
        visMeshHD.nodes = meshHD.nodes; %make copy of mesh that only contains nodes and elements for visualization purposes
        visMeshHD.elements = meshHD.elements;
end

% Visualize in coordinate space
PlotMeshSurface(visMeshHD,pM);view([0,0]) 
%view should be as if you're looking at the mesh from above and on the right side
%mesh will also be angled so that the top left corner is in the top center of the plot


%% 2f. Relax optodes onto HD Mesh
tposNew_HD=scale_cap(tposNew,1.1);
m0.nodes=meshHD.nodes;
m0.elements=meshHD.elements;
spos3=tposNew_HD(1:Ns,:);
dpos3=tposNew_HD((Ns+1):end,:);
tposNew_HD=gridspringfit_ND2(m0,rad,spos3,dpos3);

% Visualize mesh with relaxed optodes
PlotMeshSurface(visMeshHD,pM);Draw_Foci_191203(tposNew_HD,paramsFoci)
view(0,0) % Posterior view


%% Update Pad info Structure 
info.optodes.spos3=tposNew_HD(1:Ns,:);
info.optodes.dpos3=tposNew_HD((Ns+1):end,:);
m=0;
for d=1:Nd
    for s=1:Ns
        m=m+1;
        info.pairs.r3d(m)=norm(tposNew_HD(s,:)-tposNew_HD(d+Ns,:));
        info.pairs.r3d(m+(Ns*Nd))=info.pairs.r3d(m);
    end
end

% Histogram of number of measurements per SD vs radius between SD pairs
figure;histogram(info.pairs.r3d,1000);xlabel('R_S_D');ylabel('N_m_e_a_s');

% Save pad file
gridname=padname;
save(['Pad_',meshname,'_',gridname,'_', '.mat'],'info')
pM.reg = 0;


%% PREPARE! --> mesh with grid array in same file set for NIRFAST
mesh=meshHD;
mesh=PrepareMeshForNIRFAST(mesh,[meshname,'_',gridname],tposNew_HD);
PlotMeshSurface(mesh,pM);PlotSD(mesh.source.coord(1:Ns,:),...
    mesh.source.coord((Ns+1):end,:),'render',gcf);
view([0,0]) %view mesh from behind with an angled top down view

% One last visualization check...
% A portion of the mesh will be removed for this specific visualization
% This is so you can evaluate whether the optodes are placed too deeply 
% within the mesh
m3=CutMesh(mesh,intersect(find(mesh.nodes(:,3)>0),find(mesh.nodes(:,1)>0)));
[Ia,Ib]=ismember(m3.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
m3.region=mesh.region(Ib);
PlotMeshSurface(m3,pM);PlotSD(mesh.source.coord(1:Ns,:),...
    mesh.source.coord((1+Ns):end,:),'render',gcf);
view([0,0]) % View mesh from behind 


%% Section 3
%% 3a. Calculate Sensitivity Profile
% Set flags
flags.tag=[padname,'_on_',meshname];
flags.gridname=gridname;
flags.meshname=meshname;
flags.head='info';
flags.info=infoT1;                  % Your T1 info file
flags.gthresh=1e-5;                 % Voxelation threshold in G
flags.voxmm=2;                      % Voxelation resolution (mm)
flags.labels.r1='csf';              % Regions for optical properties
flags.labels.r2='white';
flags.labels.r3='gray';
flags.labels.r4='bone';
flags.labels.r5='skin';
flags.op.lambda=[750,850];          % Wavelengths (nm)
flags.op.mua_skin=[0.0170,0.0190];  % Baseline absorption
flags.op.mua_bone=[0.0116,0.0139];
flags.op.mua_csf=[0.0040,0.0040];
flags.op.mua_gray=[0.0180,0.0192];
flags.op.mua_white=[0.0167,0.0208];
flags.op.musp_skin=[0.74,0.64];     % Baseline reduced scattering coeff
flags.op.musp_bone=[0.94,0.84];
flags.op.musp_csf=[0.3,0.3];
flags.op.musp_gray=[0.8359,0.6726];
flags.op.musp_white=[1.1908,1.0107];
flags.op.n_skin=[1.4,1.4];          % Index of refration
flags.op.n_bone=[1.4,1.4];
flags.op.n_csf=[1.4,1.4];
flags.op.n_gray=[1.4,1.4];
flags.op.n_white=[1.4,1.4];
flags.srcnum=Ns;                    % Number of sources
flags.t4=eye(4);                    % T1/dim to MNI atlas *** change this to register your vol to atlas
flags.t4_target='MNI'; % string
flags.makeA=1; % don't make A, just make G
flags.Hz=0;
if flags.Hz, flags.tag = [flags.tag,'FD']; end

% Run makeAnirfast to get sensitivity matrix
Ti=tic;[A,dim,Gsd]=makeAnirfaster(mesh,flags); % size(A)= [Nwl, Nmeas, Nvox]
disp(['<makeAnirfast took ',num2str(toc(Ti))])

%Gsd vs Rsd provides a simulated light fall-off
figure;semilogy(info.pairs.r3d(1:(Ns*Nd)),Gsd,'*');
legend({num2str(flags.op.lambda')})
xlabel('R_S_D [mm]');ylabel('G_S_D');xlim([0,100])


%% 3b. Package metadata and save A
[Nwl,Nmeas,Nvox]=size(A);
A=reshape(permute(A,[2,1,3]),Nwl*Nmeas,Nvox);

% Place spatial information about light model in info.tissue structure
info.tissue.dim=dim;
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;

temp_pairs=struct;
temp_pairs.Src=info.pairs.Src;
temp_pairs.Det=info.pairs.Det;
temp_pairs.NN=info.pairs.NN;
temp_pairs.WL=info.pairs.WL;
temp_pairs.lambda=info.pairs.lambda;
temp_pairs.Mod=info.pairs.Mod;
temp_pairs.r2d=info.pairs.r2d;
temp_pairs.r3d=info.pairs.r3d;
info.pairs=temp_pairs;

save(['A_',flags.tag,'.mat'],'A','info', 'Gsd', '-v7.3') %save A


%% 3c. Update bookkeeping in A to save space and align with system data
% Get A with all SD separations within 5cm
keep=info.pairs.r3d<=50; % for sparse density pad, only save
% measurements where SD separation is within 50mm (5cm)
temp=struct;
temp.Src=info.pairs.Src(keep);
temp.Det=info.pairs.Det(keep);
temp.NN=info.pairs.NN(keep);
temp.WL=info.pairs.WL(keep);
temp.lambda=info.pairs.lambda(keep);
temp.Mod=info.pairs.Mod(keep);
temp.r2d=info.pairs.r2d(keep);
temp.r3d=info.pairs.r3d(keep);
info.pairs=temp;
A=A(keep,:); % A with all SD within 5cm

% Save A with only 1st to 5th NN
save(['A_r3d_lth_50_',flags.tag,'.mat'],'A','info', 'Gsd', '-v7.3');


%% 3d. Visualize aspects of sensitivity profile
t1=affine3d_img(mask,infoT1,dim,eye(4)); % put anatomical volume in dim space

% Visualize single source detector pair, Source 1 and Detector 24
keep=info.pairs.WL==2 & info.pairs.Src==1 & info.pairs.Det==24; % SD pair set here
sdV=squeeze(A(keep,:));              % Single meas pair
sdV=Good_Vox2vol(sdV',dim);
sdV=sdV./max(sdV(:));
sdV=log10(1e3.*sdV);               % Top 2 orders of magnitude
pA.Scale=0.8*max(sdV(:));pA.PD=1;pA.Th.P=0;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,sdV)

% Make FFR
keep=(info.pairs.WL==2 & info.pairs.r3d<=40);
a=squeeze(A(keep,:));
iA=Tikhonov_invert_Amat(a,0.01,0.1);
iA=smooth_Amat(iA,dim,5);           % 5 = smoothing parameter
ffr=makeFlatFieldRecon(a,iA);

% Save flat field reconstruction to A matrix file
save(['A_r3d_lth_50_',flags.tag,'.mat'], 'ffr', '-append');

% Visualize FFR
ffrV=Good_Vox2vol(ffr,dim);
ffrV=ffrV./max(ffrV(:));
pA.Scale = 1;
PlotSlices(t1,dim,pA,ffrV);

%% Visualize on the fsLR Surface
load([dataout, '/fs_LR/', [pt,'_', res, '_ctx.mat']])

%Visualization parameters
pA.Scale = 0.8*max(ffrV(:));pA.Th.P = 0;pA.Th.N = -pA.Th.P;pA.PD = 0;pA.view = 'post'; 

% Pial Surface
pA.ctx = 'std';
PlotInterpSurfMesh(ffrV, Anat.CtxL, Anat.CtxR, dim, pA);

% Inflated Surface
pA.ctx = 'inf';
PlotInterpSurfMesh(ffrV, Anat.CtxL, Anat.CtxR, dim, pA);

% Very Inflated Surface
pA.ctx = 'vinf';
PlotInterpSurfMesh(ffrV, Anat.CtxL, Anat.CtxR, dim, pA);

% Flat Surface
pA.view = 'lat'; pA.ctx = 'flat'; pA.reg = 1;
PlotInterpSurfMesh(ffrV, Anat.CtxL, Anat.CtxR, dim, pA);


%% Visualize alignment of LD mesh, HD mesh, array, and cortices
% Name of 'orig' cortical mesh file output from fsLR 
cd([dataout,'/fs_LR']);
fn_mesh = [pt,'_164k_ctx']; % 
% Load cortical mesh
load([fn_mesh]);

% Visualize cortical mesh inside of low-density head mesh 
pA0l=struct;
figure('Color','k','Position',[500,100,1050,1000])
pM1.fig_handle=gca;pM1.FaceColor=[0.5,0.5,0.5];pM1.FaceAlpha = 1;pM1.EdgeColor = 'none';
if isfield(pM1, 'Cmap') pM1 = rmfield(pM1, 'Cmap'); end
pM1.SpecularStrength=0.25;pM1.DiffuseStrength=0.5;pM1.BGC = [0,0,0];
PlotMeshSurface(meshHD,pM1)                     % HD mesh
pA0l.fig_handle=gca;
pA0l.FaceColor=[0.5,0.5,0.5];pA0l.EdgeColor='none';pA0l.BGC = [0,0,0];
pA0l.AmbientStrength=0.25;pA0l.DiffuseStrength=0.25;
pA0l.SpecularStrength=0.025;
PlotMeshSurface(Anat.CtxL,pA0l)                 % Cortical Surfaces
PlotMeshSurface(Anat.CtxR,pA0l)
pA2.fig_handle=gca;
pA2.FaceAlpha=0;
pA2.EdgeAlpha=0.5;
pA2.EdgeColor=[0.25,0.25,0.25];%1,1,1].*0.25;
pA2.BGC = [0,0,0];
PlotMeshSurface(meshLD,pA2)                     % LD mesh
set(gcf,'Color','k')
view([0,0]) % display in posterior view
Draw_Foci_191203(cat(1,info.optodes.spos3,info.optodes.dpos3),paramsFoci)  % Array
axis off

%Display in lateral view
view([90,0])

