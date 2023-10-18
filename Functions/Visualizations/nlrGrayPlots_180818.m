function nlrGrayPlots_180818(nlrdata,info)
%
% This function generates a gray plot figure for measurement pairs
% for just clean WL==2 data. It is assumed that the input data is
% nlrdata that has already been filtered and resampled. The data is grouped
% into info.pairs.r2d<20, 20<=info.pairs.r2d<30, and 30<=info.pairs.r2d<40.

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

[Nm,Nt]=size(nlrdata);
LineColor='w';
BkndColor='k';

if isfield(info,'GVTD'),Nrt=size(info.GVTD,1);end
% M=max(abs(nlrdata(:)));
wl=unique(info.pairs.lambda(info.pairs.WL==2));
figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.5],...
    'Color',BkndColor);

if isfield(info, 'MEAS')
    if ~isfield(info.MEAS, 'GI')
        info.MEAS.GI=ones(size(info.pairs.Src,1),1);
    end
else
    info.MEAS = struct;
    info.MEAS.GI=ones(size(info.pairs.Src,1),1);
end

%% Prepare data and imagesc together
keep.d1=info.MEAS.GI & info.pairs.r2d<20 & info.pairs.WL==2;
keep.d2=info.MEAS.GI & info.pairs.r2d>=20 & info.pairs.r2d<30 &...
            info.pairs.WL==2;
keep.d3=info.MEAS.GI & info.pairs.r2d>=30 & info.pairs.r2d<40 &...
            info.pairs.WL==2;

SepSize=round((sum(keep.d1)+sum(keep.d2)+sum(keep.d3))/50);
data1=cat(1,squeeze(nlrdata(keep.d1,:)),nan(SepSize,Nt),....*-M
        squeeze(nlrdata(keep.d2,:)),nan(SepSize,Nt),... .*-M   
        squeeze(nlrdata(keep.d3,:)));    
    
M=nanstd((data1(:))).*3;


%% Line Plot of DVARS
if isfield(info,'GVTD')
    subplot(2,1,1,'Position',[0.1,0.75,0.8,0.2])
    plot([1:Nrt],info.GVTD(:),'r','LineWidth',2);
    xlim([1,Nrt])
    set(gca,'Color',BkndColor,'XColor',LineColor,'YColor',LineColor)
    title(['GVTD'],'Color','w')
    ylabel('a.u.');
    subplot(2,1,2,'Position',[0.1,0.05,0.8,0.6])
end


%% Gray Plot data
imagesc(data1,[-1,1].*M)
hold on

% Plot synchs
DrawColoredSynchPoints(info);
Npix=size(data1,1); 
% for j=1:length(info.paradigm.synchpts)    % Draw synch pt bars
%     if ismember(j,info.paradigm.Pulse_1)
%         plot([1,1].*info.paradigm.synchpts(j),[1,Npix],'r','LineWidth',2)
%     elseif ismember(j,info.paradigm.Pulse_2)
%         plot([1,1].*info.paradigm.synchpts(j),[1,Npix],'b','LineWidth',2)
%     elseif ismember(j,info.paradigm.Pulse_3)
%         plot([1,1].*info.paradigm.synchpts(j),[1,Npix],'m','LineWidth',2)
%     elseif ismember(j,info.paradigm.Pulse_4)
%         plot([1,1].*info.paradigm.synchpts(j),[1,Npix],'g','LineWidth',2)
%     end
% end

%     pause
% Plot separators
dz1=length(keep.d1);
dz2=length(keep.d2);
dz3=length(keep.d3);
dzT=dz1+dz2+dz3+2*SepSize;
% dz=dz1+0.5;
% % rectangle('Position',[0,dz,Nt+1,SepSize],'FaceColor','k',...
% %     'EdgeColor','none')
% dz=dz1+dz2+SepSize+0.5;
% rectangle('Position',[0,dz,Nt+1,SepSize],'FaceColor','k',...
%     'EdgeColor','none')
%     pause
% Add labels
title(['\Delta',num2str(wl),' nm'],'Color',LineColor);
h1=text('String','Rsd: [1,20) mm','Units','Normalized','Position',...
    [-0.04,(dzT-0.45*dz1)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');
h2=text('String','Rsd: [20,30) mm','Units','Normalized','Position',...
    [-0.04,(dz3+SepSize+0.6*dz2)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');
h3=text('String','Rsd: [30,40) mm','Units','Normalized','Position',...
    [-0.04,(0.60*dz3)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');

set(gca,'XTick',[],'YTick',[],'Box','on','Color','w');

set(gcf,'Color',BkndColor)
colormap(gray(1000))