clear 
close all
clc

%% load t1w, t2w, pdw, bi-mask
[T1,infoT1] = LoadVolumetricData(['nihpd_asym_05-08_t1w'], [],'nii');
[T2,infoT2] = LoadVolumetricData(['nihpd_asym_05-08_t2w'], [],'nii');
[T3,infoT3] = LoadVolumetricData(['nihpd_asym_05-08_pdw'], [],'nii');
[T4,infoT4] = LoadVolumetricData(['nihpd_asym_05-08_mask'], [],'nii');
T1 = T1./max(T1(:));
T2 = T2./max(T2(:));
T3 = T3./max(T3(:));
T4 = T4./max(T4(:));
p.slices_type='coord';p.slices=[-0.5,0.5,0.5];p.CH=0;
p.Cmap='jet';p.Scale=1;p.Th.P=0;p.Th.N=-p.Th.P;p.PD=1;p.BG=[0,0,0];
PlotSlices(T1,infoT1,p) % T1
PlotSlices(T2,infoT2,p) % T2
PlotSlices(T3,infoT3,p) % PDW
PlotSlices(T4,infoT4,p) % Mask

%% extract head
head_T1_bw=ExtractHead(T1,0.17);
head_T1_bw=imfill(head_T1_bw,'holes');
head=+(head_T1_bw);
p.Cmap='jet';p.Scale=1;p.Th.P=0;p.Th.N=-p.Th.P;p.PD=1;p.BG=[0,0,0];
PlotSlices(head,infoT1,p) 

%% Get mask including CSF, GM, WM
brain2 = T2.*T4;
brain2(brain2<0.5) = 0;  % fix the errors in mask
brain2 = imfill(brain2);
PlotSlices(brain2,infoT1,p) % T2 and mask
brain3 = T3.*T4.*+(brain2>0);
PlotSlices(brain3,infoT1,p) % PDW and mask
brain1 = T1.*T4.*+(brain2>0);
PlotSlices(brain1,infoT1,p) % T1 and mask

%% extract WM
CSF1 = +(brain1<0.4 & brain1>0);
PlotSlices(CSF1,infoT1,p)

%%
brainfilter =  +(brain3>0.9);
PlotSlices(brainfilter,infoT1,p)

%%
WM0 = +((brain1>0.7 & brain1<0.8) | (brain2>0.85 & brain2<0.89))-CSF1;
PlotSlices(WM0,infoT1,p)

%%
WM = WM0.*brainfilter;
PlotSlices(WM,infoT1,p)

%% extract CSF in brain
CSFin = brain1-WM;
CSFin(CSFin>0.6) = 0;
CSFin = +(CSFin>0)++((T4-brain1)==1);
PlotSlices(CSFin,infoT1,p)

%% extract GM
GM = +((brain1-CSFin-WM)>0);
PlotSlices(GM,infoT1,p)

%% extract surface including skin, skull, CSF
%% extract CSF
head(:,:,1:20) = 0;
surf1 = T1.*head-T4;
surf2 = T2.*head-T4;
surf3 = T3.*head-T4;
PlotSlices(surf1,infoT1,p)
PlotSlices(surf2,infoT1,p)
PlotSlices(surf3,infoT1,p)
%%
CSFout = +(surf2>0.85);
PlotSlices(CSFout,infoT1,p)

%% extract skin and skull
extra = +((surf2-CSFout)>0);
extra = imfill(extra);
PlotSlices(extra,infoT1,p)

%%
CSFout = +((surf2-extra)>0);
PlotSlices(CSFout,infoT1,p)
CSF = CSFin+CSFout;

%% make 4-layer head model
mask=zeros(size(T1));
mask(extra==1)=4;
mask(GM==1)=3;
mask(WM==1)=2;
mask(CSF==1)=1;
p.Cmap='jet';p.Scale=4;p.Th.P=0;p.Th.N=-p.Th.P;p.PD=1;p.BG=[0,0,0];
PlotSlices(mask,infoT1,p)
load('affine_nihpd_atlas');
affine_0508 = affine_nihpd_atlas;
mask = affine3d_img(mask,infoT1,infoT1,affine_0508,'nearest');

save('6mo_4layer_Mask.mat','mask','infoT1')