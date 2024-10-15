function FreeSurfer_Recon_all_PSHM(subjid,fn,dataroot, mode)
%
% This function runs FreeSurfer's recon-all scripts to extract the cortical
% ribbon from a native space t1w MRI.
% Inputs:
%   subjid   subject ID used on MRI data
%   fn      full filename (w/o extensions) for native-space t1 to be sent
%               through recon-all. Typically this is the cleanest possible
%               t1w volume - it has the best registration metric of those
%               t1w present in the atlas directory.
%   dataroot    the full path to the folder tree containing the anatomical
%               data. Typically, this is a path that ends with the folder
%               '/atlas/'.

%% Prepare Inputs
if ~exist('mode','var')
    mode = '4dfp';
end
%% Prepare T1 for FreeSurfer
disp('<<< Preparing T1 for FreeSurfer')% Transform to nifti
if strcmp(mode, '4dfp')
    [status,response]=system(['nifti_4dfp -n ',fn,' ',fn,''])
end
cd /usr/local/freesurfer/subjects/  % Go to FreeSurfer director 
[status,response]=system(['recon-all -i ',...% Organize files for a new subject
    dataroot,fn,'".nii" -subjid ',subjid,''])

%% Run recon-all
disp('<<< Running FreeSurfer recon-all')
t=tic;                              
[status,response]=system(['recon-all -subjid ',subjid,' -all -clean'])
toc(t)

%% Run freesurfer_to_fs_LR
cd /usr/local/freesurfer/subjects/ 
% t1=tic;
% evalc(['system([''/usr/local/freesurfer/subjects/fsLR_ES/',...
%     'freesurfer_to_fs_LR.sh ',subjid,'''])'])
% disp('<<< Pipeline to Caret fs_LR')
% toc(t1)

%% Move data to dataroot tree
[status,response]=system(['mkdir ',dataroot,'freesurfer/']);
[status,response]=system(['mv ',subjid,'/* ',dataroot,'freesurfer/']);
evalc(['system([''mkdir ',dataroot,'fsLR/''])']);
% evalc(['system([''mv fs_LR_output_directory/',subjid,'/* ',...
%     dataroot,'fsLR/''])']);
