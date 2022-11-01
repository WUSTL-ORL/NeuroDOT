function Plot_RawData_Cap_levels(data,info,params)
%
% This function plots various views of light levels and noise levels.
% The layout of the figure includes all-source-by-all-detector grids on the
% top left and top center showing the temporal standard deviation in Y. 
% These plots help illuminate if there is something wrong with a given
% source or detector in terms of the electronics.
% On the top right is a light fall-off plot. On the bottom left is a cap
% layout of the mean across neighbors of the mean light levels 
% (mean in time, AKA Phi_0) over a given radius for each source and 
% detector on the cap. On the bottom center is a similar figure but 
% showing mean across neighbors of the standard deviation in Y. These 2
% plots help point out which optode positions are good and which are bad, 
% which can often be due to poor coupling at the scalp. On the bottom left
% is a 'Good measurements' plot of green lines for measurements passing the
% noise threshold (default is std(Y)<0.075). The neighbors of focus for
% these green lines can be the same as those for the light levels and noise
% levels (default), or they can be set to be different.  e.g., if the user
% wants Cap-based mean light and noise levels for source-detector pairs
% between 10-20 mm but good measurements between 20-30mm, set 
% params.rlimits=[10,20]; and params.GMrlimits=[20,30];
% Because typically, the second WL is the longer wavelength, those
% measurements have higher SNR and are used as the default for the
% cap-based visualizations. Passing params.WL can be used to change which 
% WL is of focus.
%
%% Parameters and Initialization
warning off
if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end
if ~isfield(params,'LLmax'),params.LLmax=-1;end  % log10 for cap light levels cMap
if ~isfield(params,'LLmin'),params.LLmin=-3;end % log10
dr=params.LLmax-params.LLmin;
if ~isfield(params,'WL'),params.WL=2;end
if ~isfield(params,'rlimits')
    params.rlimits=[1,20];
end
Rlimits=params.rlimits;
if ~isfield(params,'GMrlimits'),params.GMrlimits=Rlimits;end
if ~isfield(params, 'mode'),params.mode = 'good';end
if ~isreal(data),data=abs(data);end

if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    %     params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    params.fig_handle=figure('Units','Normalized',...
        'Position',[0.1,0.1,0.75,0.8],'Color','k');
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end

Nm = sum(info.pairs.WL==params.WL);
Ns = length(unique(info.pairs.Src));
Nd = length(unique(info.pairs.Det));
cs = unique(info.pairs.WL); % WLs.

params.PD = 1;
params.Th.P = 0;
params.DR = 1000;
params.Cmap.P = 'hot';
params.dimension = '2D';
cmThreshMult=1.5;


dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);

wls=unique(info.pairs.lambda);
Nwls=length(wls);
leg=cell(Nwls,1);
for j=1:Nwls, leg{j}=[num2str(wls(j)),' nm'];end

keep=(info.pairs.lambda==min(wls));
dList=repmat([1:Nd],Ns,1);
measFull=[repmat([1:Ns]',Nd,1),dList(:)];
[Ia,Ib]=ismember(measFull,[info.pairs.Src(keep),info.pairs.Det(keep)],'rows');
Ib(Ib==0)=[];
sdFull=nan(Ns,Nd);

NNrIgnore=min(info.pairs.NN(info.pairs.r2d>35));

cMap=jet(1000);

%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Initial math
Phi_0 = mean(abs(data), 2);
Y=logmean(abs(data));
stdY=std(Y,[],2);
stdYcf=stdY;
stdYcf(info.pairs.NN>=NNrIgnore)=0;

%% std(Y) WL 1
if isfield(info,'paradigm')
    if isfield(info.paradigm, 'synchpts')
        NsynchPts = length(info.paradigm.synchpts); % set timing of data
        if NsynchPts > 2
            tF = info.paradigm.synchpts(end);
            t0 = info.paradigm.synchpts(2);
        elseif NsynchPts == 2
            tF = info.paradigm.synchpts(2);
            t0 = info.paradigm.synchpts(1);
        else
            tF = size(data, 2);
            t0 = 1;
        end
        stdY=std(Y(:, t0:tF),[],2);
    else
        stdY=std(Y,[],2);
    end
else
    stdY=std(Y,[],2);
end
sdFull(Ia)=stdY(Ib);
subplot(2,3,1,'Units','Normalized','Position',[0.025,0.55,0.3,0.4]) 
imagesc(sdFull,[0,0.2]);
colormap(gca,cMap);
title(['\sigma (Y) ',num2str(wls(1)),' nm'],'Color','w');xlabel('Detector')
ylabel('Source');set(gca,'XColor','w','YColor','w');axis square;


%% std(Y) WL 2
sdFull(Ia)=stdY(Ib+Nm);
subplot(2,3,2,'Position',[0.30,0.55,0.3,0.4])
imagesc(sdFull,[0,0.2]);
p0=get(gca,'Position');
title(['\sigma (Y) ',num2str(wls(2)),' nm'],'Color','w');xlabel('Detector')
colormap(gca,cMap);
cb=colorbar('Ticks',[0,0.075,0.1,0.2],...
    'TickLabels',{'0','0.075','\sigma(Y)','0.2'},...
    'Color','w','Position',[[p0(1)+p0(3)+0.005,p0(2)+0.075,0.005,0.2]]);
set(gca,'XColor','w','YColor','w','YTickLabel','');axis square;


%% Light level fall off
subplot(2,3,3,'Position',[0.70,0.55,0.27,0.4])
params.fig_handle=gca;
Plot_LightLevel_FallOff(data,info,params);
xM=min([100,max(info.pairs.r2d)+5]);
xlim([0,xM]);


%% Light levels at optode
subplot(2,3,4,'Position',[0.025,0.025,0.3,0.4])
params.fig_handle=gca;

for s = 1:Ns
    keep = (info.pairs.Src == s)  &  info.pairs.r2d>Rlimits(1,1) &...
        info.pairs.r2d<=Rlimits(1,2) & info.pairs.WL==params.WL;
    Srcval(s) = log10(mean(Phi_0(keep)));
end

for d = 1:Nd
    keep = (info.pairs.Det == d)  &  info.pairs.r2d>Rlimits(1,1) &...
        info.pairs.r2d<=Rlimits(1,2) & info.pairs.WL==params.WL;
    Detval(d) = log10(mean(Phi_0(keep)));
end

% Absolute color mapping
Srcval = Srcval - (params.LLmax - dr);
Detval = Detval - (params.LLmax - dr);
params.Scale=dr;
params.BG=[0,0,0];
[SDRGB, CMAP] = applycmap([Srcval(:); Detval(:)], [], params);
SrcRGB = SDRGB(1:Ns, :);
DetRGB = SDRGB(Ns+1:end, :);
params.mode = 'textpatch';
PlotCapData(SrcRGB, DetRGB, info, params);
title(['\mu (\Phi_0); ',num2str(Rlimits(1,1)),'<R_s_d<',...
    num2str(Rlimits(1,2)),'; ',...
    num2str(max(info.pairs.lambda)),' nm'],'Color','w')
colormap(gca,'hot');caxis([0,1])
cMin=params.LLmax-dr;
cbMLL=colorbar('Location','southoutside','Color','w',...
    'LimitsMode','manual','Limits',[0,1],...
    'Ticks',[0,0.5,1],...,(params.LLmax+cMin)/2,params.LLmax],...
    'TickLabels',{[num2str(10^cMin)],'a.u.',...
    [num2str(10^params.LLmax)]});


%% Noise levels at optode
subplot(2,3,5,'Position',[0.35,0.025,0.3,0.4])
params.fig_handle=gca;

for s = 1:Ns
    keep = (info.pairs.Src == s)  &  info.pairs.r2d>=Rlimits(1,1) &...
        info.pairs.r2d<=Rlimits(1,2) & info.pairs.WL==params.WL;
    Srcval(s) = (mean(stdY(keep)));
end

for d = 1:Nd
    keep = (info.pairs.Det == d)  &  info.pairs.r2d>=Rlimits(1,1) &...
        info.pairs.r2d<=Rlimits(1,2) & info.pairs.WL==params.WL;
    Detval(d) = (mean(stdY(keep)));
end

% Absolute color mapping
params.Scale=params.bthresh*cmThreshMult;
[SDRGB, CMAP] = applycmap([Srcval(:); Detval(:)], [], params);
SrcRGB = SDRGB(1:Ns, :);
DetRGB = SDRGB(Ns+1:end, :);

params.LineColor=[0,0,0]; % to see bad S/D identity
PlotCapData(SrcRGB, DetRGB, info, params);
title(['\mu (\sigma); ',num2str(Rlimits(1,1)),'<R_s_d<',...
    num2str(Rlimits(1,2)),'; ',...
    num2str(max(info.pairs.lambda)),' nm'],'Color','w')
colormap(gca,'hot');caxis([0,1])
cbMNL=colorbar('Location','southoutside','Color','w',...
    'LimitsMode','manual','Limits',[0,1],...[0,params.bthresh],...
    'Ticks',[0,0.5,cmThreshMult^(-1),1],...params.bthresh/2,params.bthresh],...
    'TickLabels',{'0','th %',num2str(params.bthresh*100),...
    num2str(params.bthresh*cmThreshMult*100)});


%% Good Measurements
subplot(2,3,6,'Position',[0.68,0.025,0.275,0.35])
info.MEAS.STD = std(Y,[],2);
info.MEAS.GI = info.MEAS.STD<=params.bthresh;
Pgm.fig_handle=gca;
Pgm.rlimits=params.GMrlimits;
PlotCapGoodMeas(info, Pgm);