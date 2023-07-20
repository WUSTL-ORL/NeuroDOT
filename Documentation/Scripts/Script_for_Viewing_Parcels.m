%% Script for Viewing Parcels on Surfaces

load('IM_Gordon_2014_333_Parcels.mat','IM','Parcels','Parcel_Nets')
load('Conte69_on_TT_32k.mat')


%% View all parcels on Conte69
% (1) Add parcel nodes to Cortices
Anat.CtxL.data=Parcels.CtxL;
Anat.CtxR.data=Parcels.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=jet(333);
params.TC=1;
params.ctx='inf';           % 'std','inf','vinf'
params.view='dorsal';       % 'dorsal','post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);


%% View Parcels with coloring based on network
Anat.CtxL.data=Parcel_Nets.CtxL;
Anat.CtxR.data=Parcel_Nets.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=IM.cMap;
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'
params.view='lat';        % also, 'post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);


%% View only parcels within a specific network
net='None';
keep=ismember(IM.Nets,net);
cMap=IM.cMap;
cMap(~keep,:)=0.5;
cMap(keep,:)=[1,0,0];
params.Cmap.P=cMap;
params.TC=1;
params.ctx='std';         % also, 'std','inf','vinf'
params.view='dorsal';        % also, 'post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
