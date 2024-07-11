function info = Plot_RawData_Metrics_II_DQC(data,info,params)
%
% This function generates a single-page report that includes:
%   zoomed raw time trace
%   light fall off as a function of Rsd
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
warning off
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



%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(logmean(data), info, params.bthresh,params);
end


%% Make Figure
params.fig_handle=figure('Units','Normalized',...
    'Position',[0.05,0.05,0.6,0.8],'Color','k');


%% Light level fall off
Phi_0=mean(data,2);
M=ceil(max(log10(data(:))));
yM=10^M;
subplot(3,2,[2,4],'Position',[0.55,0.42,0.4,0.55]) 
r=info.pairs.r3d(info.pairs.lambda==wls(1));
semilogy(r,reshape(Phi_0,Nm,[]),'.');
axis([0,100,1e-6,yM])
xlabel('Source-Detector Separation ( mm )','Color','w')
ylabel('{\Phi_0} ( {\mu}W )','Color','w')
set(gca,'XColor','w','YColor','w','Xgrid','on','Ygrid','on','Color','k')
legend(leg,'Color','w')

if isfield(info.MEAS,'Clipped')
    hold on;
    keep=info.MEAS.Clipped==1;
    semilogy(info.pairs.r3d(keep),Phi_0(keep),'xw');
    legend(cat(1,leg,'Clipped'),'Color','k','TextColor','w')
end

m=ceil(min(log10(Phi_0(:))))-2;
Phi_0_to_plot=reshape(Phi_0,Nm,[]);
Keep_range = find(r<20);

info.DQ_metrics.WL1.min_val = min(Phi_0_to_plot(Keep_range,1));
info.DQ_metrics.WL1.min_log = log10(info.DQ_metrics.WL1.min_val);
info.DQ_metrics.WL1.max_val = max(Phi_0_to_plot(Keep_range,1));
info.DQ_metrics.WL1.max_log = log10(info.DQ_metrics.WL1.max_val);
info.DQ_metrics.WL1.range = info.DQ_metrics.WL1.max_log - info.DQ_metrics.WL1.min_log;
info.DQ_metrics.WL1.WL = wls(1);

info.DQ_metrics.WL2.min_val= min(Phi_0_to_plot(Keep_range,2));
info.DQ_metrics.WL2.min_log = log10(info.DQ_metrics.WL2.min_val);
info.DQ_metrics.WL2.max_val = max(Phi_0_to_plot(Keep_range,2));
info.DQ_metrics.WL2.max_log = log10(info.DQ_metrics.WL2.max_val);
info.DQ_metrics.WL2.range = info.DQ_metrics.WL2.max_log - info.DQ_metrics.WL2.min_log;
info.DQ_metrics.WL2.WL = wls(2);


%% <fft> for max wavelength at 2 distances
lmdata = logmean(data);

subplot(3,2,5,'Position',[0.05,0.05,0.4,0.3]) 
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


%% Noise Histogram
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
keep=(info.pairs.lambda==min(wls));
dList=repmat([1:Nd],Ns,1);
measFull=[repmat([1:Ns]',Nd,1),dList(:)];
[Ia,Ib]=ismember(measFull,[info.pairs.Src(keep),info.pairs.Det(keep)],'rows');
Ib(Ib==0)=[];
keep=intersect(find(info.pairs.r3d>=min(params.rlimits(:)) & ...
        info.pairs.r3d<=max(params.rlimits(:))),Ib);
subplot(3,2,6,'Position',[0.55,0.05,0.4,0.3])    
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
subplot(3,2,[1,3],'Position',[0.05,0.42,0.4,0.55])
if sum(keep)
semilogy([ti:tf]./fr,squeeze(data(keep,ti:tf)))
end
set(gca,'XColor','w','YColor','w','Color','k')
ttmin = 1*10^(floor(log10(mean(Phi_0)))); % Adjust y-axis limits based on input data
ttmax = 1*10^(floor(log10(mean(Phi_0)))+1);
axis([[ti,tf]./fr,ttmin,ttmax])
title(['\Phi(t) ',num2str(wls(2)),' nm, GI: Rsd from ',...
    num2str(min(params.rlimits(:))),' - ',...
    num2str(max(params.rlimits(:))),' mm'],'Color','w')
xlabel('Time [sec]')

%