function GrayPlot_TimeTrace_Rsd_by_WL(data,info,params)
%
% This function plots time traces of data that has already been log-meaned
% using plot with different rows corresponding to different .
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
if ~isfield(params,'rlimits'),params.rlimits=[1,20;21,30;31,40];end
if ~isfield(params, 'yscale'),params.yscale = 'linear';end
lambdas=unique(info.pairs.lambda,'stable');
WLs=unique(info.pairs.WL,'stable');

params.fig_handle=figure('Position',[100 100 1050 780],'Color','w');

%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(data, info, params.bthresh);
end


%% Draw data
for j=1:Nwl
    for k=1:3
        subplot(3,Nwl,((k-1)*Nwl)+j)
        
        keep=info.pairs.r2d>=params.rlimits(k,1) & ...
            info.pairs.r2d<=params.rlimits(k,2) & ...
            info.pairs.WL==j & info.MEAS.GI;
         
        m=max(max(abs(data(keep,:))));
        imagesc(data(keep,:),[-1,1].*m);
        set(gca,'XLimSpec','tight'), xlabel('Time (samples)'),
        ylabel([{['Meas Idx']};{['within group']}]);
        colorbar%('Location','northoutside');
        colormap gray
        if isfield(info.pairs,'lambda')
            title(['Good ',num2str(lambdas(j)),' nm, Rsd:',...
                num2str(params.rlimits(k,1)),'-',...
                num2str(params.rlimits(k,2)),' mm'])
        else
            title(['Good WL ## ',num2str(WLs(j)),' nm, Rsd:',...
                num2str(params.rlimits(k,1)),'-',...
                num2str(params.rlimits(k,2)),' mm'])
        end
        
        if isfield(info,'paradigm') % Add in experimental paradigm timing
            DrawColoredSynchPoints(info,0);
        end
    end
end