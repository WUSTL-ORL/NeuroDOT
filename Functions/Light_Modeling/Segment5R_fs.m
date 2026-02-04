function [mask,params]=Segment5R_fs(T1n,T2,bm,aseg,info,params)

% This function takes as input information from FreeSurfer's volumetric
% segmentation and then finishes the job with the outer head sections.
% T1n:      B-field normalized T1
% T2:       t2w coaligned with the T1n in mpr1 space
% bm:       brainmask from FS T1
% aseg:     segmented volume from FS
% asegStat: data from ascii file about regions in aseg
%
%
% FS does a great job in brain, not so much outside of it.
% CSF:  T1 very dark, T2 very bright
% WM:   T1 bright, T2 dark
% GM:   T1 dark, T2 bright
% Scalp: T1 brighter, T2 darker
% Skull: (compact) T1 dark, T2 dark

%% Set parameters
if ~exist('info','var'),info.acq='transverse';end
if ~exist('params','var'),params=struct;end
if ~isfield(params,'Rsa'),params.Rsa=2;end % skull dilate: skull thickness
if ~isfield(params,'Rba2'),params.Rba2=2;end % bknd dilate(2): final scalp boundary
if ~isfield(params,'Np'),params.Np=10;end               % buffer size
if ~isfield(params,'Head_th'),params.Head_th=0.1;end    % Head thresh T1
if ~isfield(params,'Head_th2'),params.Head_th2=0.1;end    % Head thresh T2
if ~isfield(params,'head_bg_seg_type'),params.head_bg_seg_type='T1T2';end % Head bkgnd seg type
if ~isfield(params,'Skull_th'),params.Skull_th=0.25;end % Skull thresh
if ~isfield(params,'T2_bias_kernel'),params.T2_bias_kernel=8; end % csf T1 thresh
if ~isfield(params,'CSF_T1_th'),params.CSF_T1_th=0.175; end % csf T1 thresh
if ~isfield(params,'CSF_T2_th'),params.CSF_T2_th=0.2750; end % csf T2 thresh
if ~isfield(params,'killY'),params.killY=[1:20]; end    % ventral removal
if ~isfield(params,'smooth_gauss'),params.smooth_gauss=1.2; end    % Gaussian smoothing

SE_sa = strel('ball', params.Rsa, 0); % skull dilate: skull thickness
SE_sa2 = strel('ball', params.Rba2, 0);% bknd dilate(2): final scalp boundary

[Nx,Ny,Nz]=size(T1n);


%% zero pad to be safe
disp('<<< Adding buffer to volumes')
T1n=Zero_Pad_Vol(T1n,info,params.Np);
T2=Zero_Pad_Vol(T2,info,params.Np);
bm=Zero_Pad_Vol(bm,info,params.Np);
aseg=Zero_Pad_Vol(aseg,info,params.Np);


%% Normalize volumes
disp('<<< Filling holes in head')
T1n=T1n./max(T1n(:));
T2=T2./max(T2(:));
bm(bm>0)=1;
bm=imfill(bm,'holes');


%% Separate head from background
disp('<<< Extracting Head')
switch params.head_bg_seg_type
    case 'T1T2'
        head_T1_bw=ExtractHead(T1n,params.Head_th);  % Extract T1 head shape
        head_T2_bw=ExtractHead(T2,params.Head_th2);  % Extract T2 head shape
        head_T1_bw=imfill(head_T1_bw,'holes');
        head_T2_bw=imfill(head_T2_bw,'holes');
        head=+(head_T1_bw.*head_T2_bw);
    case 'T1'
        head_T1_bw=ExtractHead(T1n,params.Head_th);  % Extract T1 head shape
        head_T1_bw=imfill(head_T1_bw,'holes');
        head=+(head_T1_bw);
    case 'T2'
        head_T2_bw=ExtractHead(T2,params.Head_th2);  % Extract T2 head shape
        head_T2_bw=imfill(head_T2_bw,'holes');
        head=+(head_T2_bw);
end
head=imfill(head,'holes');
head=+(smooth3(head,'gaussian',[7,7,7],2)>0.5);    % Smooth head
bknd=1-head;
head_T1_data=T1n.*head;                     % masked T1


%% Get fs brain segmentation
mask1=SegmentFSAseg(aseg);
GM_bw=mask1.GM.*bm;
WM_bw=mask1.WM.*bm;
CSF_bw=mask1.CSF.*bm;
clear mask1;


%% Clean up CSF and GM in brain mask
leftover=+((bm-GM_bw-WM_bw-CSF_bw)>0);
csf2=(((T2)>=params.CSF_T2_th).*((T1n)<=params.CSF_T1_th)).*leftover;
gm2=(T1n>params.CSF_T1_th).*leftover;
CSF_bw=CSF_bw+csf2;
GM_bw=+((GM_bw+gm2)>0);
to_kill=((T2<params.CSF_T2_th).*(T1n<=params.CSF_T1_th)).*leftover;
bm(to_kill==1)=0;


%% Brain mask separation from Scalp+Skull
disp('<<< Correcting Tissue Masks')
SS_bw=1-(bm+bknd);
SS_T1=SS_bw.*head_T1_data;
SS_T2=SS_bw.*T2;


%% Scalp Skull separation
Skull_bw=(SS_T2<params.Skull_th).*(SS_T1<params.Skull_th).*SS_bw;
Skull_bw=imdilate(Skull_bw,SE_sa);                  % Dilate Skull
Skull_bw=SingleContiguousRegion(Skull_bw);          % Skull SCR
Skull_bw=imfill(Skull_bw,'holes');                  % Fill holes in skull
Skull_bw=imerode(Skull_bw,SE_sa);                   % Erode Skull
Skull_bw=SingleContiguousRegion(Skull_bw).*SS_bw;   % Skull SCR
Scalp_bw=SS_bw-Skull_bw;                            % Scalp


%% Set values to region numbers: [1,2,3,4,5]=[csf,wm,gm,sk,sc]
disp('<<< Smoothing head and setting mask values')

head_bknd=1-head;
head_bknd=+(imdilate(head_bknd,SE_sa2)>0);
head_bknd=+(imgaussfilt3(head_bknd,params.smooth_gauss)>0.5);    % Smooth head
head_T1_bw=1-head_bknd;

mask=head.*5;
mask(Skull_bw==1)=4;
mask(Scalp_bw==1)=5;
mask(CSF_bw==1)=1;
mask(WM_bw==1)=2;
mask(GM_bw==1)=3;

% mask(:,1:65,:)=0;
mask=mask.*head_T1_bw;
mask=mask(:,:,((params.Np+1):(end-params.Np)));% fix z
mask=mask(:,((params.Np+1):(end-params.Np)),:);% fix y
mask=mask(((params.Np+1):(end-params.Np)),:,:);% fix x

% Remove extra neck stuff
switch info.acq(1)
    case 't'
mask(:,:,params.killY)=0;
    case 's'
mask(:,params.killY,:)=0;
end
