function mask=Segment5R_fs_noT2(T1n,bm,aseg,info,params)

% This function takes as input information from FreeSurfer's volumetric
% segmentation and then finishes the job with the outer head sections.
% T1n:      B-field normalized T1
% T2:       t2w coaligned with the T1n in mpr1 space
% bm:       brainmask from FS
% aseg:     segmented volume from FS
% asegStat: data from ascii file about regions in aseg


%% Set parameters
if ~exist('params','var'),params=struct;end
if ~isfield(params,'Rbr'),params.Rbr=5;end
if ~isfield(params,'Rsa'),params.Rsa=2;end
if ~isfield(params,'Rba'),params.Rba=3;end
if ~isfield(params,'Rba2'),params.Rba2=2;end
if ~isfield(params,'Np'),params.Np=10;end               % buffer size
if ~isfield(params,'Head_th'),params.Head_th=0.1;end    % Head thresh
if ~isfield(params,'Skull_th'),params.Skull_th=0.25;end % Skull thresh
if ~isfield(params,'CSF_T1_th'),params.CSF_T1_th=0.250; end % csf thresh
if ~isfield(params,'killY'),params.killY=[1:20]; end    % ventral removal

SE_br = strel('ball', params.Rbr, 0); % brain dilate: outer skull from brain
SE_sa = strel('ball', params.Rsa, 0); % skull dilate: skull thickness
SE_ba = strel('ball', params.Rba, 0); % bknd dilate: skull thickness
SE_sa2 = strel('ball', params.Rba2, 0);% bknd dilate: final scalp boundary
% [Nx,Ny,Nz]=size(T1n);


%% Add buffer
disp('<<< Adding buffer to volumes')
T1n=Zero_Pad_Vol(T1n,params.Np);
bm=Zero_Pad_Vol(bm,params.Np);
aseg=Zero_Pad_Vol(aseg,params.Np);


%% Normalize volumes
disp('<<< Filling holes in head')
T1n=T1n./max(T1n(:));
bm(bm>0)=1;
bm=imfill(bm,'holes');


%% Separate head from background
disp('<<< Extracting Head')
head_T1_bw=ExtractHead(T1n,params.Head_th);   % Extract T1 head shape
head_T1_bw=imfill(head_T1_bw,'holes');
head=+(head_T1_bw); %.*head_T2_bw);           % mask FOV based on T2 FOV
head=imfill(head,'holes');
head=+(imgaussfilt3(head,2)>0.5);             % Smooth head
bknd=1-head;
head_T1_data=T1n.*head;                     % masked T1



%% Get fs brain segmentation
mask1=SegmentFSAseg(aseg);
GM_bw=mask1.GM;
WM_bw=mask1.WM;
CSF_bw=mask1.CSF;
clear mask1;


%% Clean up CSF and GM in Brain mask
leftover=+((bm-GM_bw-WM_bw-CSF_bw)>0);
csf2=((T1n)<=params.CSF_T1_th).*leftover;
gm2=(T1n>params.CSF_T1_th).*leftover;
CSF_bw=CSF_bw+csf2;
GM_bw=+((GM_bw+gm2)>0);
to_kill=((T1n<=params.CSF_T1_th)).*leftover;
bm(to_kill==1)=0;


%% Brain mask
disp('<<< Correcting Tissue Masks')
SS_bw=1-(bm+bknd);
SS_T1=SS_bw.*head_T1_data;


%% Scalp Skull separation
Skull_bw=(SS_T1<params.Skull_th).*SS_bw;
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
head_bknd=+(imgaussfilt3(head_bknd,1.2)>0.5);             % Smooth head
head_T1_bw=1-head_bknd;

mask=head.*5;
% mask=zeros(size(GM_bw));
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
