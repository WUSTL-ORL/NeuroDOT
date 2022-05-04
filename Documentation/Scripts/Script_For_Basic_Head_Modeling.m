%% Script for model generation given an initial segmented vokumetric mask and optode set

% Volumetric segmentation is currently not supported, but there are many
% options available (e.g., FreeSurfer). The starting point to this script
% is a segmented volume and a Pad file, both of which are contained within
% the Support_Files folder of the toolbox. To generate your own head model
% and light model, please just adapt the steps outlined below to your own
% volumetric mask and source-detector grid design.


%% Load a Segmented Volume and generate head mesh

[mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');
p.Cmap='jet';p.Scale=5;p.Th.P=0;p.Th.N=-p.Th.P;p.PD=1;p.BG=[0,0,0];
PlotSlices(mask,infoT1,p)       % Visualize the segmented mask: 
        % Note, PlotSlices is an interactive plot. 
        % Hit q or the middle mouse button to quit


% Parameters for generating your mesh
meshname=['Example_Mesh'];      % Provide a name for your mesh name here
param.facet_distance=0.75;   % Node position error tolerance at boundary
param.facet_size =0.75;      % boundary element size parameter
param.cell_size=1.4;        % Volume element size parameter
param.info=infoT1;
param.Offset=[0,0,0];
param.r0=5;                  % nodes outside of mask must be set to scalp==5;
param.CheckMeshQuality=0;
param.Mode=0;
tic;mesh=NirfastMesh_Region(mask,meshname,param);toc

pM.orientation='coord';pM.Cmap.P='gray';
PlotMeshSurface(mesh,pM)
%%%% IF you get an error, this is due to NIRFAST using an old version of
%%%% mode. Go into NIRFAST/toolbox, change mode.m to modeOLD.m. Then the
%%%% correct version of mode will be used and this visualization will work.

% Put coordinates back in true space
mesh.nodes=change_space_coords(mesh.nodes,infoT1,'coord');
PlotMeshSurface(mesh,pM)

% Look at the inside
m1.nodes=mesh.nodes;m1.elements=mesh.elements;
m2=CutMesh(m1,find(m1.nodes(:,2)>0));
[Ia,Ib]=ismember(m2.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
m2.region=mesh.region(Ib);
PlotMeshSurface(m2,pM);


%% Load initial Pad file depending on data structure:
padname='AdultV24x28';
load(['Pad_',padname,'.mat']); % Pad file
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3);
Ns=size(info.optodes.spos3,1);
Nd=size(info.optodes.dpos3,1);
rad=info;
PlotSD(tpos(1:Ns,:),tpos((Ns+1):end,:),'norm');
PlotCap(info)


%% View cap position with and without mesh
if ~exist('mesh','var')
    [mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');
    mesh=load_mesh(meshname);
    mesh.nodes=change_space_coords(mesh.nodes,infoT1,'coord');
    PlotMeshSurface(mesh,pM)
end

PlotMeshSurface(mesh,pM);PlotSD(tpos(1:Ns,:),tpos((Ns+1):end,:),'norm',gcf);


%% Move grid from arbitrary location to approximate target on mesh

% Adjust parameters by hand, update position, check
tpos2=tpos;
dx=0;
dy=-110;
dz=-55;
dS=1.1;
dxTh=-90;
dyTh=0;
dzTh=0;

tpos2=rotate_cap(tpos2,[dxTh,dyTh,dzTh]);
tpos2(:,1)=tpos2(:,1)+dx;
tpos2(:,2)=tpos2(:,2)+dy;
tpos2(:,3)=tpos2(:,3)+dz;
tpos2=scale_cap(tpos2,dS);

PlotMeshSurface(mesh,pM);PlotSD(tpos2(1:Ns,:),tpos2((Ns+1):end,:),'norm',gcf);


%% Relax grid on head
m0.nodes=mesh.nodes;
m0.elements=mesh.elements;
spos3=tpos2(1:Ns,:);
dpos3=tpos2((Ns+1):end,:);
tposNew=gridspringfit_ND2(m0,rad,spos3,dpos3);

info.optodes.spos3=tposNew(1:Ns,:);
info.optodes.dpos3=tposNew((Ns+1):end,:);
info.tissue.infoT1=infoT1;
info.tissue.affine=eye(4);
info.tissue.affine_target='MNI';
save(['Pad_',padname,'_on_',meshname,'.mat'],'info','tposNew'); % Pad file


%% View optode positions
PlotMeshSurface(mesh,pM);PlotSD(tposNew(1:Ns,:),tposNew((Ns+1):end,:),'render',gcf);

%% PREPARE! --> mesh and grid in coord space
gridname=padname;
mesh=PrepareMeshForNIRFAST(mesh,[meshname,'_',gridname],tposNew);
PlotMeshSurface(mesh,pM);PlotSD(mesh.source.coord(1:Ns,:),...
    mesh.source.coord((Ns+1):end,:),'render',gcf);

% One last visualization check...
m3=CutMesh(mesh,intersect(find(mesh.nodes(:,3)>0),find(mesh.nodes(:,1)>0)));
[Ia,Ib]=ismember(m3.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
m3.region=mesh.region(Ib);
PlotMeshSurface(m3,pM);PlotSD(mesh.source.coord(1:Ns,:),...
    mesh.source.coord((1+Ns):end,:),'render',gcf);



%% Calculate Sensitivity Profile
% Parameters
flags.tag=[padname,'_on_',meshname,'_test'];
flags.gridname=gridname;
flags.meshname=meshname;
flags.head='info';
flags.info=infoT1;                  % Your T1 info file
flags.gthresh=1e-3;                 % Voxelation threshold in G
flags.voxmm=2;                      % Voxelation resolution (mm)
flags.labels.r1='csf';             % Regions for optical properties
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

Ti=tic;[A,dim]=makeAnirfast(mesh,flags); % size(A)= [Nwl, Nmeas, Nvox]
disp(['<makeAnirfast took ',num2str(toc(Ti))])


%% Package data and save A with only 1st through 5th nearest neighbors to save space

nn12345=1; % for 96x92 cap

[Nwl,Nmeas,Nvox]=size(A);
A=reshape(permute(A,[2,1,3]),Nwl*Nmeas,Nvox);

info.tissue.dim=dim;
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;

if nn12345
keep=info.pairs.NN<=5;
temp=table;
temp.Src=info.pairs.Src(keep);
temp.Det=info.pairs.Det(keep);
temp.NN=info.pairs.NN(keep);
temp.WL=info.pairs.WL(keep);
temp.lambda=info.pairs.lambda(keep);
temp.Mod=info.pairs.Mod(keep);
temp.r2d=info.pairs.r2d(keep);
temp.r3d=info.pairs.r3d(keep);
info.pairs=temp;
A=A(keep,:);
end

save(['A_',flags.tag,'.mat'],'A','info','-v7.3')



%% Visualize aspects of sensitivity profile
t1=affine3d_img(mask,infoT1,dim,eye(4)); % put anatomical volume in dim space

keep=info.pairs.WL==2 & info.pairs.Src==1 & info.pairs.Det==6;
foo=squeeze(A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo,dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)

% FFR
keep=(info.pairs.WL==2 & info.pairs.r2d<=40);
a=squeeze(A(keep,:));
iA=Tikhonov_invert_Amat(a,0.01,0.1);
iA=smooth_Amat(iA,dim,3);
ffr=makeFlatFieldRecon(a,iA);

fooV=Good_Vox2vol(ffr,dim);
fooV=fooV./max(fooV(:));
pA.PD=1;pA.Scale=1;pA.Th.P=1e-2;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)
