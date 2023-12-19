%% OVERVIEW
% 1. start with Mask and Pad
% 2. make a Low Density mesh of the mask
% 3. use AlignMe to place the Pad on the LD mesh
% 4. remove unneeded components of the mask
% 5. make a HD mesh of the scooped mask
% 6. run AlignMe on the LD-aligned Pad and the HD mesh
% 7. update the info file and name it specifically for this light model
% 8. prepare mesh for nirfast
% 9. set flags and run makeAnirfast
% 10. package and save A with info etc
% 11. view single measurement, Gsd, and flat field reconstruction
% 12. view LD mesh, HD mesh, array, and cortex together


%% NOTE
% When running this script, please DO NOT HIT THE GREEN TRIANGLE "RUN" BUTTON
% Please run individual code sections by highlighting and evalutating the highlighted selection
% Or, you can hit "ctrl + enter" on each section of code once the section's background turns beige


%% Pathing and Parameters
%Set Output Directory
%Change the path on the following line to whatever output directory you want to use
outputdir = ''; % put path to directory in between quotes
cd(outputdir)

%Parameter Initialization
p.Cmap='jet'; p.Scale=5; p.Th.P=0; p.Th.N=-p.Th.P; p.PD=1; p.BG=[0,0,0];
pM.orientation='coord'; pM.Cmap.P='gray'; pM.EdgeColor='none';
pS = pM; pS = rmfield(pS, 'orientation'); % make copy of params, to be used with PlotSlices


%% Load a Segmented Volume and generate head mesh
[mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111_nifti'],[],'nii');
PlotSlices(mask,infoT1,p)       % Visualize the segmented mask: 
        % Note, PlotSlices is an interactive plot. 
        % Hit q or the middle mouse button to quit


%% Generate low density head mesh
% Parameters for generating your mesh
meshname=['LD_Mesh'];      % Provide a name for your mesh name here
param.facet_distance=5;    % Node position error tolerance at boundary, changed to 5 from .75 12/2/19
param.facet_size=3;        % boundary element size parameter
param.cell_size=5;         % Volume element size parameter
param.info=infoT1;
param.Offset=[0,0,0];
param.r0=5;                % nodes outside of mask must be set to scalp==5;
param.CheckMeshQuality=0;
param.Mode=1;              % make simple mesh with no region labels
tic;meshLD=NirfastMesh_Region(mask,meshname,param);toc

PlotMeshSurface(meshLD,pM) %visualize in index (non-negative) space
view([135,30]) %view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is

%%%% IF you get an error, this is due to NIRFAST using an old version of
%%%% mode. Go into NIRFAST/toolbox, change mode.m to modeOLD.m. Then the
%%%% correct version of mode will be used and this visualization will work.

% Put coordinates back in true space
meshLD.nodes=change_space_coords(meshLD.nodes,infoT1,'coord');
PlotMeshSurface(meshLD,pM) %visualize in coordinate space
view([135,30])%view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is
%now the mesh should be centered at (0,0,0)

% Look at the inside
m1.nodes=meshLD.nodes;m1.elements=meshLD.elements;
m2=CutMesh(m1,find(m1.nodes(:,2)>0)); %cut LD mesh in half
[Ia,Ib]=ismember(m2.nodes,meshLD.nodes,'rows');Ib(Ib==0)=[];
PlotMeshSurface(m2,pM); %visualize back half of LD mesh
view([135,30])%view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is
%now the mesh should be centered at (0,0,0)


%% Load initial Pad file depending on data structure:
padname='AdultV24x28';
load(['Pad_',padname,'.mat']); % Pad file
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3);
Ns=size(info.optodes.spos3,1);
Nd=size(info.optodes.dpos3,1);
rad=info;

paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1));
paramsFoci.color(1,:) = [1 0.4 0.6]; % pink for s1
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for d1
paramsFoci.radius = 2; %set radius of spheres to 2

%Visualize Pad
Draw_Foci_191203(tpos, paramsFoci);view([-40,30]) %3D plot of pad
PlotCap(info) %2D plot of pad

% Visualize LD mesh with optode positions 
% Optodes have not been relaxed onto mesh yet
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tpos, paramsFoci);
view([0,90])%view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is
%you should be able to see a few optodes peeking out of the mesh
%the array has yet to be lined up correctly, which you'll do in the next section with hummingbird


%% AlignMe Section:  Move grid from arbitrary location to approximate target on mesh, Relax grid on head and view
% Location of atlas fiducials                    Mesh2EEG output
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
myapp=AlignMe_2020b(ds);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew = ds.dO.tpos2_relaxed;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew,paramsFoci)
view(0,0) %posterior view

% Get affine matrix that can be used to transform FROM participant space TO MNI space
if isfield(ds.dO, 'affineTform') %if mesh scaled, affineTform field will exist, save it to workspace
    affine_Subj2MNI = [ds.dO.affineTform, zeros(3,1)];
    save('affine_matrix_Subject_to_MNI.mat', 'affine_Subj2MNI')
else %otherwise, set affine_Subj2MNI to eye(4)
    affine_Subj2MNI = eye(4);
end


%% Generate smaller mask based on optode locations
maskCrop=IceCreamScoopCap4Mask(infoT1,tposNew,mask,40); %40 = mm of mesh to scoop out
PlotSlices(maskCrop,infoT1,pS) %visualize smaller mask


%% Generate High Density Head Mesh (Can take up to 15 min)
%If you get an error when running NirfastMesh_Region
%try chagning the mesh name and clearing your output directory of all files
meshname=['HD_Mesh0'];      % Provide a name for your mesh name here, if making multiple meshes, provide a different name for each mesh
param.facet_distance=2.0;   % Node position error tolerance at boundary
param.facet_size=0.8;       % boundary element size parameter
param.cell_size=1.0;        % Volume element size parameter
param.Mode=0;               % Set Mode=0 to make nirfast compliant mesh
param.info=infoT1;          % Make sure info is infoT1 like with LD mesh
param.Offset=[0,0,0];
tic;meshHD=NirfastMesh_Region(maskCrop,meshname,param);toc

%make copy of mesh that only contains nodes and elements for visualization purposes
visMeshHD.nodes = meshHD.nodes;
visMeshHD.elements = meshHD.elements;
PlotMeshSurface(visMeshHD,pM);view([70,60]) %Visualize in index space
%view should be as if you're looking at the mesh from above and on the right side
%mesh will also be angled so that the top left corner is in the top center of the plot

% Put coordinates back in coordinate space
meshHD.nodes=change_space_coords(meshHD.nodes,infoT1,'coord'); %for actual mesh
visMeshHD.nodes=change_space_coords(visMeshHD.nodes,infoT1,'coord'); %for mesh that's visualized
PlotMeshSurface(visMeshHD,pM);view([70,60]) %Visualize in coordinate space
%view should be as if you're looking at the mesh from above and on the right side
%mesh will also be angled so that the top left corner is in the top center of the plot


%% Relax optodes onto HD Mesh
% Create an instance of  DataStorage class
ds_HD = DataStorage();

% Input structure
ds_HD.dI.tpos = tposNew;       
ds_HD.dI.mesh = meshHD;
ds_HD.dI.pad = info;
ds_HD.dI.pM = pM;
ds_HD.dI.paramsFoci = paramsFoci;
ds_HD.dI.Ns = Ns;
ds_HD.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds_HD.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds_HD.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application,
% passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
myapp = AlignMe_2020b(ds_HD);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew_HD = ds_HD.dO.tpos2_relaxed;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(visMeshHD,pM);Draw_Foci_191203(tposNew_HD,paramsFoci)
view(0,0) %posterior view


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

% histogram of number of measurements per SD vs radius between SD pairs
figure;histogram(info.pairs.r3d,1000);xlabel('R_S_D');ylabel('N_m_e_a_s');

% save pad file
gridname=padname;
save(['Pad_',meshname,'_',gridname,'_', '.mat'],'info')


%% PREPARE! --> mesh with grid array in same file set for NIRFAST
mesh=meshHD;
mesh=PrepareMeshForNIRFAST(mesh,[meshname,'_',gridname],tposNew_HD);
PlotMeshSurface(mesh,pM);Draw_Foci_191203(tposNew_HD, paramsFoci);
view([0,0]) %posterior view


% One last visualization check...
% A portion of the mesh will be removed for this specific visualization
% This is so you can evaluate whether the optodes are placed too deeply within the mesh
m3=CutMesh(mesh,intersect(find(mesh.nodes(:,3)>0),find(mesh.nodes(:,1)>0)));
[Ia,Ib]=ismember(m3.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
m3.region=mesh.region(Ib);
PlotMeshSurface(m3,pM);Draw_Foci_191203(tposNew_HD, paramsFoci);
view([135,40]) %view mesh from behind with an angled top-down view



%% Calculate Sensitivity Profile
% Set flags
flags.tag=[padname,'_on',meshname,'_test'];
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
flags.op.n_skin=[1.4,1.4];          % Index of refraction
flags.op.n_bone=[1.4,1.4];
flags.op.n_csf=[1.4,1.4];
flags.op.n_gray=[1.4,1.4];
flags.op.n_white=[1.4,1.4];
flags.srcnum=Ns;                    % Number of sources
flags.t4=affine_Subj2MNI;           % Affine matrix for going from subject-specific space to MNI space
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


%% Package data and save A
[Nwl,Nmeas,Nvox]=size(A);
A=reshape(permute(A,[2,1,3]),Nwl*Nmeas,Nvox);

% Place spatial information about light model in info.tissue structure
info.tissue.dim=dim;
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;

save(['A_',flags.tag,'.mat'],'A','info','-v7.3') %save A


%% Visualize aspects of sensitivity profile
t1=affine3d_img(mask,infoT1,dim,eye(4)); % put anatomical volume in dim space

% Visualize Single SD pair (S1 D6)
keep=info.pairs.WL==2 & info.pairs.Src==1 & info.pairs.Det==6; % SD pair set here
foo=squeeze(A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)

% Visualize Single SD pair (S2 D2)
keep=info.pairs.WL==2 & info.pairs.Src==2 & info.pairs.Det==2; % SD pair set here
foo=squeeze(A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)

% FFR
keep=(info.pairs.WL==2 & info.pairs.r2d<=40);
a=squeeze(A(keep,:));
iA=Tikhonov_invert_Amat(a,0.01,0.1);
iA=smooth_Amat(iA,dim,5); %5 = smoothing parameter
ffr=makeFlatFieldRecon(a,iA);

fooV=Good_Vox2vol(ffr,dim);
fooV=fooV./max(fooV(:));
pA.PD=1;pA.Scale=1;pA.Th.P=1e-2;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)


%% Visualize alignment of LD mesh, HD mesh, array, and cortices
%  View alignment
Cortical_mesh = load(['MNI164k_big.mat']);
Anat.CtxL = Cortical_mesh.MNIl;
Anat.CtxR = Cortical_mesh.MNIr;
pA0l=struct;
figure('Color','k','Position',[500,100,1050,1000])
pM1.fig_handle=gca;pM1.FaceColor=[1,1,1].*0.5;pM1.EdgeColor = 'none';
pM1.SpecularStrength=0.25;pM1.DiffuseStrength=0.5;
% PlotMeshSurface(meshHD,pM1)                        % HD mesh
pA0l.fig_handle=gca;
pA0l.FaceColor=[1,1,1].*0.5;pA0l.EdgeColor='none';
pA0l.AmbientStrength=0.25;pA0l.DiffuseStrength=0.25;
pA0l.SpecularStrength=0.025; pA0l.ctx = 'inf';
PlotMeshSurface(Anat.CtxL,pA0l)                 % Cortical Surfaces
PlotMeshSurface(Anat.CtxR,pA0l)
pA2.fig_handle=gca;
pA2.FaceAlpha=0;
pA2.EdgeAlpha=0.5;
pA2.EdgeColor=[1,1,1].*0.25;
PlotMeshSurface(meshLD,pA2)                     % LD mesh
set(gca,'Color','k')
view([90,0]) % display in Lateral view
Draw_Foci_191203(cat(1,info.optodes.spos3,info.optodes.dpos3),paramsFoci)  % Array
axis off

% For posterior view, uncomment the following line
% view([0,0])

