function Plot_Measurement_PER(data,info,StimF)
%
% This function calculates the Phase-Encoded-Retinotopy analyses and then
% creates a figure of source-detector rectangles showing the total variance,
% Magnitude, Phase, and Significant magnitude at the stimulus frequency
% (StimF) for all source-detector pairs.
%
%
%% Parameters and Initialization
if isfield(info, 'MEAS')
    if istable(info.MEAS)
        info.MEAS = table2struct(info.MEAS, 'ToScalar', true);
    end
end
if isfield(info, 'pairs')
    if istable(info.pairs)
        info.pairs = table2struct(info.pairs, 'ToScalar', true);
    end
end

fs=info.system.framerate;
Trange=info.paradigm.synchpts(2):info.paradigm.synchpts(end);
Ns=max(info.pairs.Src);
Nd=max(info.pairs.Det);
lambdas=unique(info.pairs.lambda);
wls=unique(info.pairs.WL);
Nwl=length(wls);
WL=2;
keep=(info.pairs.WL==min(wls));
Nm=sum(keep);
adjSD=Ns*Nd==size(info.pairs.Src);
cMap=cat(1,[0,0,0],jet(1000));
cMapP=cat(1,[0,0,0],hot(1000));
cMapPh=cat(1,[0,0,0],hsv(1000));

params.fig_handle=figure('Color','k','Units','Normalized',...
    'Position',[0.15,0.25,0.8,0.25]);


%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(data, info, 0.075);
end


%% Check for non-full measurement list
if adjSD
    dList=repmat([1:Nd],Ns,1);
    measFull=[repmat([1:Ns]',Nd,1),dList(:)];
    [Ia,Ib]=ismember(measFull,[info.pairs.Src(keep),info.pairs.Det(keep)],'rows');
    Ib(Ib==0)=[];
    Ib=Ib+Nm;
else
    Ia=1:Nm;
    Ib=find(info.pairs.WL==WL);
end
sdFull=nan(Ns,Nd);


%% Calculate PER stuff
[Mag0,Phase0,SigMag0]=Phase_Enc_Ret(data(:,Trange),fs,StimF);


%% Variance
subplot(1,4,1);
Vsd=var(data,[],2);
sdFull(Ia)=Vsd(Ib).*info.MEAS.GI(Ib);
imagesc(sdFull);
title(['\sigma^2 (Y) ',num2str(lambdas(WL)),' nm'],'Color','w');
colormap(gca,cMap);
cb=colorbar('Color','w');
xlabel('Detector');ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;


%% Magnitude at stim freq (relative scale)
subplot(1,4,2);
sdFull(Ia)=Mag0(Ib).*info.MEAS.GI(Ib);
sdFull=sdFull./max(sdFull(:));
imagesc(sdFull,[0,1]);
title(['Mags = |FFT (Y; ',num2str(StimF),' Hz)|'],'Color','w');
colormap(gca,cMap);
cb=colorbar('Ticks',[0,0.5,1],'TickLabels',{'0','50%','100%'},'Color','w');
xlabel('Detector');ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;


%% Significant Magnitudes at stim freq (relative scale)
subplot(1,4,3);
sdFull(Ia)=SigMag0(Ib).*info.MEAS.GI(Ib); % SigMag are p-values
imagesc(sdFull,[-log10(0.05),10]);
title(['Significance of Mags (-log_1_0(p))'],'Color','w');
colormap(gca,cMapP);
cb=colorbar('Color','w');
xlabel('Detector');ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;


%% Phase at stim freq 
subplot(1,4,4);
sdFull(Ia)=Phase0(Ib).*info.MEAS.GI(Ib); % SigMag are p-values
sdFull(sdFull==0)=nan;
imagesc(sdFull,[-pi,pi]);
title(['Phase at ',num2str(StimF),' Hz'],'Color','w');
colormap(gca,cMapPh);
cb=colorbar('Ticks',[-pi,0,pi],'TickLabels',{'-\pi','0','\pi'},'Color','w');
xlabel('Detector');ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;