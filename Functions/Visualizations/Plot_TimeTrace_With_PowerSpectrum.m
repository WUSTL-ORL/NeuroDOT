function Plot_TimeTrace_With_PowerSpectrum(data,info,params)
%
% This function plots time traces of data that has already been log-meaned
% using plot (top row), imagesc (middle row), and the power spectrum of the
% meas across kept measurements (bottom row).
%



%% Parameters and Initialization
Nwl=length(unique(info.pairs.WL));
[Nm,Nt]=size(data);
if isfield(info.system,'init_framerate')
    fr=info.system.init_framerate;
else
    fr=info.system.framerate;
end
    
dt=1/fr;                          
t=[1:1:Nt].*dt;

if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end
if ~isfield(params,'rlimits'),params.rlimits=[1,40];end
if ~isfield(params, 'yscale'),params.yscale = 'log';end
lambdas=unique(info.pairs.lambda,'stable');
WLs=unique(info.pairs.WL,'stable');

params.fig_handle=figure('Position',[100 100 1050 780],'Color','w');

%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(data, info, params.bthresh);
end


%% Draw data
for j=1:Nwl
    keep=info.pairs.r2d>=params.rlimits(1,1) & ...
        info.pairs.r2d<=params.rlimits(1,2) & ...
        info.pairs.WL==j & info.MEAS.GI;
    
    subplot(3,Nwl,j)
    plot(data(keep,:)'); 
    m=max(max(abs(data(keep,:))));
    set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
    ylabel('log(\phi/\phi_0)');
    mData=squeeze(mean(data(keep,:),1));
    hold on
    plot(mData,'k','LineWidth',2)
    ylim([[-1,1].*m])
    if isfield(info.pairs,'lambda')
        title(['Good ',num2str(lambdas(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'])
    else
        title(['Good WL ## ',num2str(WLs(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'])
    end    
    
if isfield(info,'paradigm') % Add in experimental paradigm timing
    DrawColoredSynchPoints(info,0);
end
    
    subplot(3,Nwl,j+Nwl) 
    imagesc(data(keep,:),[-1,1].*m.*0.5); 
    colorbar('Location','northoutside');
    colormap gray
    xlabel('Time (samples)');ylabel('Measurement #')
    if isfield(info.pairs,'lambda')
        title(['Good ',num2str(lambdas(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'])
    else
        title(['Good WL ## ',num2str(WLs(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'])
    end    
    
if isfield(info,'paradigm') % Add in experimental paradigm timing
    DrawColoredSynchPoints(info,0);
end

    
    subplot(3,Nwl,j+2*Nwl);
    [ftmag0,ftdomain] = fft_tts(data(keep,:),...
        info.system.framerate); % Generate average spectrum
    ftmag=rms(ftmag0,1);
    semilogx(ftdomain,ftmag);
    xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])

    
end