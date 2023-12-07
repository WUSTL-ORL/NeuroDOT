function DQ_metrics = Plot_RawData_Metrics_I_DQC(data,info,params)
%
% This function generates a single-page report that includes:
%   light fall off as a function of Rsd
%   SNR plots vs. mean light level
%   source-detector mean light-level plots
%   Power spectra for 830nm at 2 Rsd
%   Histogram for measurement noise
% 
%
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

%% Parameters and Initialization
if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end
if ~isfield(params,'rlimits')
    params.rlimits=[1,20;21,35;36,43];
elseif size(params.rlimits,1)==1
    params.rlimits=cat(1,params.rlimits,[21,35;36,43]);
elseif size(params.rlimits,1)==2
    params.rlimits=cat(1,params.rlimits,[36,43]);
end
if isfield(info.system,'init_framerate')
    fr=info.system.init_framerate;
else
    fr=info.system.framerate;
end
if ~isfield(params,'logfft'),params.logfft=0;end
if ~isfield(params,'LFO_GI'),params.LFO_GI=0;end

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

wls=unique(info.pairs.lambda);
Nwls=length(wls);
leg=cell(Nwls,1);
for j=1:Nwls, leg{j}=[num2str(wls(j)),' nm'];end

Ns=max(info.pairs.Src);
Nd=max(info.pairs.Det);
[Nm,Nt]=size(data);
if ~isreal(data),data=abs(data);end
Nm=Nm/Nwls;

if Nt<(60/fr)
    ti=1;
    tf=Nt;
elseif Nt
    ti=round(Nt/2)-round(5*fr);
    if ti<1, ti=1;end
    tf=round(Nt/2)+round(5*fr);
    if tf>=Nt, tf=Nt;end
end

% Add in more exact correction for APD Sensitivity (~6e6 V/W) to convert Mags
% to Voltage to OptPower(rms). Also include estimate for theoretical noise
% floor via NEP*(sqrt(bandwidth))
if ~isfield(params,'Input_Refer_Power'),params.Input_Refer_Power=1;end
if params.Input_Refer_Power
    if isfield(params,'Sens')
        APDsens=params.Sens;
    else
        APDsens=6e6; % ~6 MV/W is sensitivity of our APDs
    end
    data=(data./APDsens).*1e6; % Input refer data to optical power in micro Watts
end

if ~isfield(params,'NEPth'),params.NEPth=0;end
if params.NEPth
    % Theoretical noise floor should look like NEP*sqrt(bw)
    NEPth=params.NEPth;
end


%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(logmean(data), info, params.bthresh,params);
end


%% Make Figure
params.fig_handle=figure('Units','Normalized',...
    'Position',[0.05,0.05,0.8,0.8],'Color','k');


%% Light level fall off
if params.LFO_GI==1
    keep=params.LFO_GI;
    keep=(sum(reshape(keep,Nm,[]),2)>=1);
else
    keep=ones(Nm,1)==1;
end
Phi_0=mean(data,2);
M=ceil(max(log10(Phi_0(:))));
yM=10^M;
m=ceil(min(log10(Phi_0(:))))-2;
ym=10^m;
subplot(3,6,[1,2,7,8],'Position',[0.05,0.42,0.28,0.55]) 
Phi_0_to_plot=reshape(Phi_0,Nm,[]);
r=info.pairs.r3d(info.pairs.lambda==wls(1));
semilogy(r(keep),Phi_0_to_plot(keep,:),'.');
axis([0,100,ym,yM])
xlabel('Source-Detector Separation ( mm )','Color','w')
ylabel('{\Phi_0 _R_M_S} ( {\mu}W )','Color','w')
set(gca,'XColor','w','YColor','w','Xgrid','on','Ygrid','on','Color','k')
legend(leg,'Color','w')

if isfield(info, 'MEAS')
    if istable(info.MEAS)
        info.MEAS = table2struct(info.MEAS, 'ToScalar', true);
    end
end

if isfield(info.MEAS, 'Clipped')
    hold on;
    if params.LFO_GI~=0
        keep = info.MEAS.Clipped & params.LFO_GI;
    else
        keep = info.MEAS.Clipped;
    end
    semilogy(info.pairs.r3d(keep), Phi_0(keep),'xw');
    leg = cat(1,leg,'Clipped');
    legend(leg, 'Color', 'k', 'TextColor', 'w');
end

if params.NEPth
    hold on
    semilogy([0:100],ones(101,1).*NEPth,'--w');
    legend(cat(1,leg,'Noise Floor Estimate'),'Color','k','TextColor','w')
end

Phi_0_to_plot=reshape(Phi_0,Nm,[]);
Keep_range = find(r<20);

DQ_metrics = struct;
DQ_metrics.min_val_WL1 = min(Phi_0_to_plot(Keep_range,1));
DQ_metrics.min_log_WL1 = log10(DQ_metrics.min_val_WL1);
DQ_metrics.max_val_WL1 = max(Phi_0_to_plot(Keep_range,1));
DQ_metrics.max_log_WL1 = log10(DQ_metrics.max_val_WL1);
DQ_metrics.range_WL1 = DQ_metrics.max_log_WL1 - DQ_metrics.min_log_WL1;

DQ_metrics.min_val_WL2= min(Phi_0_to_plot(Keep_range,2));
DQ_metrics.min_log_WL2 = log10(DQ_metrics.min_val_WL2);
DQ_metrics.max_val_WL2 = max(Phi_0_to_plot(Keep_range,2));
DQ_metrics.max_log_WL2 = log10(DQ_metrics.max_val_WL2);
DQ_metrics.range_WL2 = DQ_metrics.max_log_WL2 - DQ_metrics.min_log_WL2;

%% <fft> for max wavelength at 2 distances
lmdata = logmean(data);

subplot(3,6,13:14,'Position',[0.05,0.05,0.28,0.3])
keep=(info.pairs.lambda==max(wls)) & info.MEAS.GI & ...
    info.pairs.r3d>=params.rlimits(1,1) & ...
    info.pairs.r3d<=params.rlimits(1,2);
r1=mean(info.pairs.r3d(keep));

% WL2reg_a=mean(squeeze(lmdata(keep,:)),1);
WL2reg_a=lmdata(keep,:);
[ftmag0, ftdomain] = fft_tts(WL2reg_a,fr);
ftmag=rms(ftmag0,1);
if params.logfft
    loglog(ftdomain,ftmag,'--r','LineWidth',1);hold on
else
    semilogx(ftdomain,ftmag,'--r','LineWidth',1);hold on
end

keep=(info.pairs.lambda==max(wls)) & info.MEAS.GI & ...
    info.pairs.r3d>=params.rlimits(2,1) & ...
    info.pairs.r3d<=params.rlimits(2,2);
r2=mean(info.pairs.r3d(keep));
% WL2reg_b=mean(squeeze(lmdata(keep,:)),1);
WL2reg_b=lmdata(keep,:);
[ftmag0, ftdomain] = fft_tts(WL2reg_b,fr);
ftmag=rms(ftmag0,1);

if params.logfft
    loglog(ftdomain,ftmag,'-m','LineWidth',1);
else
    semilogx(ftdomain,ftmag,'-m','LineWidth',1);
end
xlim([1e-3,fr/2])
xlabel('Frequency [Hz]');
ylabel('|P1 [au]|');
set(gca,'XColor','w','YColor','w','Xgrid','on','Ygrid','on','Color','k')
legend([{['\mu(',num2str(max(wls)),' nm, ~',...
    num2str(r1),' mm)']};...
    {['\mu(',num2str(max(wls)),' nm, ~',...
    num2str(r2),' mm)']}],...
    'Color','w')


%% Phi_0 WL 1
keep=(info.pairs.lambda==min(wls));
dList=repmat([1:Nd],Ns,1);
measFull=[repmat([1:Ns]',Nd,1),dList(:)];
[Ia,Ib]=ismember(measFull,[info.pairs.Src(keep),info.pairs.Det(keep)],'rows');
Ib(Ib==0)=[];
sdFull=zeros(Ns,Nd);
sdFull(Ia)=Phi_0(Ib);
subplot(3,6,3,'Position',[0.37,0.66,0.14,0.35])       
imagesc(log10(sdFull),[-5,-1]);colormap(jet(1000))
title(['\Phi_0 ',num2str(wls(1)),' nm'],'Color','w');xlabel('Detector')
ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;
cMap=cat(1,[0,0,0],jet(1000));colormap(gca,cMap);


%% Phi_0 WL 2
sdFull(Ia)=Phi_0(Ib+Nm);
subplot(3,6,4,'Position',[0.52,0.66,0.14,0.35])       
imagesc(log10(sdFull),[-5,-1]);colormap(gca,cMap);
p0=get(gca,'Position');
title(['\Phi_0 ',num2str(wls(2)),' nm'],'Color','w');xlabel('Detector')
cb=colorbar('Ticks',[-5,-3,-1],'TickLabels',{'10^-^5','\muW','10^-^1'},...
    'Color','w','Position',[[p0(1)+p0(3)+0.005,p0(2)+0.075,0.005,0.2]]);
set(gca,'XColor','w','YColor','w','YTickLabel','')
axis square;colormap(gca,cMap);


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
        stdY=std(lmdata(:, t0:tF),[],2);
    else
        stdY=std(lmdata,[],2);
    end
else
    stdY=std(lmdata,[],2);
end
sdFull(Ia)=stdY(Ib);
subplot(3,6,9,'Position',[0.37,0.36,0.14,0.35])
imagesc(sdFull,[0,0.2]);
colormap(gca,cMap);
title(['\sigma (Y) ',num2str(wls(1)),' nm'],'Color','w');xlabel('Detector')
ylabel('Source');set(gca,'XColor','w','YColor','w');axis square;

%% std(Y) WL 2
sdFull(Ia)=stdY(Ib+Nm);
subplot(3,6,10,'Position',[0.52,0.36,0.14,0.35])       
imagesc(sdFull,[0,0.2]);
p0=get(gca,'Position');
title(['\sigma (Y) ',num2str(wls(2)),' nm'],'Color','w');xlabel('Detector')
colormap(gca,cMap);
cb=colorbar('Ticks',[0,0.075,0.1,0.2],...
    'TickLabels',{'0','0.075','\sigma(Y)','0.2'},...
    'Color','w','Position',[[p0(1)+p0(3)+0.005,p0(2)+0.075,0.005,0.2]]);
set(gca,'XColor','w','YColor','w','YTickLabel','');axis square;


%% SNR WL 1
snrN=nan(size(sdFull));
snrD=nan(size(sdFull));
snrN(Ia)=Phi_0(Ib);
snrD(Ia)=std(data(Ib,:),[],2);
subplot(3,6,15,'Position',[0.37,0.05,0.14,0.35])      
imagesc(log10(snrN./snrD),[0.5,2]);
colormap(gca,cMap);
title(['SNR ',num2str(wls(1)),' nm'],'Color','w');xlabel('Detector')
ylabel('Source');set(gca,'XColor','w','YColor','w');
axis square;


%% SNR WL 2
snrN(Ia)=Phi_0(Ib+Nm);
snrD(Ia)=std(data(Ib+Nm,:),[],2);
subplot(3,6,16,'Position',[0.52,0.05,0.14,0.35]);        
imagesc(log10(snrN./snrD),[0.5,2]);
colormap(gca,cMap);
p0=get(gca,'Position');
title(['SNR ',num2str(wls(2)),' nm'],'Color','w');axis square;xlabel('Detector');
cb=colorbar('Ticks',[0.5,1.25,2],'TickLabels',{'1.2','SNR','10^2'},...
    'Color','w','Position',[[p0(1)+p0(3)+0.005,p0(2)+0.075,0.005,0.2]]);
set(gca,'XColor','w','YColor','w','YTickLabel','')


%% Noise Histogram
keep=intersect(find(info.pairs.r3d>=min(params.rlimits(:)) & ...
        info.pairs.r3d<=max(params.rlimits(:))),Ib);
subplot(3,6,17:18,'Position',[0.73,0.05,0.24,0.23])    
[h1,x1]=hist(stdY(keep).*100,[0:0.5:100]);
[h2,x2]=hist(stdY(keep+Nm).*100,[0:0.5:100]);
b1=bar(x1,h1,'b');
hold on;
b2=bar(x2,h2,'g');
yl=get(gca,'Ylim');axis([0,30,0,ceil(max([h1,h2]))])
plot(ones(1,2).*params.bthresh.*100,[0,yl(end)],'r','LineWidth',2)
legend(cat(1,leg,'Threshold'),'Color','w')
xlabel('\sigma(y) [ % ]');ylabel('Measurements')
set(gca,'XColor','w','YColor','w','Color','k')
title(['Rsd from ',num2str(min(params.rlimits(:))),' - ',...
    num2str(max(params.rlimits(:))),' mm'],'Color','w')


%% Time trace bit
keep=info.pairs.r3d>=min(params.rlimits(:)) & ...
        info.pairs.r3d<=max(params.rlimits(:)) & ...
        info.MEAS.GI & info.pairs.lambda==max(wls);
subplot(3,6,[5,6,11,12],'Position',[0.72,0.36,0.25,0.605])
if sum(keep)
semilogy([ti:tf]./fr,squeeze(data(keep,ti:tf)))
end
set(gca,'XColor','w','YColor','w','Color','k')
axis([[ti,tf]./fr,1e-2,1e-1])
title(['\Phi(t) ',num2str(wls(2)),' nm, GI: Rsd from ',...
    num2str(min(params.rlimits(:))),' - ',...
    num2str(max(params.rlimits(:))),' mm'],'Color','w')
xlabel('Time [sec]')

%