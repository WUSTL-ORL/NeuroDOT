%% Load inputs
load('AlignMe_Tutorial_Inputs.mat')


%% Visualize Mesh by itself
PlotMeshSurface(mesh,pM) %visualize in coordinate space
view([135,30])%view the mesh where the nose of the head is pointing down and to the right
%this will give you a good view of where the center of the mesh is


%% Visualize optode array (pad) by itself
figure
Draw_Foci_191203(tpos, paramsFoci);view([-40,30]) %3D plot of pad
PlotCap(info) %2D plot of pad


%% Visualize mesh and array together
% Optodes have not been relaxed onto mesh yet
PlotMeshSurface(mesh,pM);Draw_Foci_191203(tpos, paramsFoci);
view([0,90])


%% Use AlignMe:  Move grid from arbitrary location to approximate target on mesh, Relax grid on head and view               
% Create an instance of our custom DataStorage HANDLE class to store variables
ds = DataStorage(); 

% Input structure
ds.dI.tpos = tpos;       
ds.dI.mesh = mesh;
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


