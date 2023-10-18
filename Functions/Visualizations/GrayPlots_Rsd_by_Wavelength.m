function GrayPlots_Rsd_by_Wavelength(data,info)
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

Nwl=length(unique(info.pairs.WL));
[Nm,Nt]=size(data);
LineColor='w';
BkndColor='k';
% Nrt=size(info.GVTD,1);
% M=max(abs(nlrdata(:)));
wl=unique(info.pairs.lambda);
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
for j=1:Nwl
    subplot(1,Nwl,j)

keep.d1=info.MEAS.GI & info.pairs.r2d<20 & info.pairs.WL==j;
keep.d2=info.MEAS.GI & info.pairs.r2d>=20 & info.pairs.r2d<30 &...
            info.pairs.WL==j;
keep.d3=info.MEAS.GI & info.pairs.r2d>=30 & info.pairs.r2d<40 &...
            info.pairs.WL==j;

SepSize=round((sum(keep.d1)+sum(keep.d2)+sum(keep.d3))/50);
data1=cat(1,squeeze(data(keep.d1,:)),nan(SepSize,Nt),....*-M
        squeeze(data(keep.d2,:)),nan(SepSize,Nt),... .*-M   
        squeeze(data(keep.d3,:)));    
    
M=nanstd((data1(:))).*3;

imagesc(data1,[-1,1].*M)
hold on

% Plot synchs
DrawColoredSynchPoints(info);

% Plot separators
dz1=length(keep.d1);
dz2=length(keep.d2);
dz3=length(keep.d3);
dzT=dz1+dz2+dz3+2*SepSize;


title(['\Delta',num2str(wl(j)),' nm'],'Color',LineColor);

if j==1
h1=text('String','Rsd: [1,20) mm','Units','Normalized','Position',...
    [-0.04,(dzT-0.45*dz1)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');
h2=text('String','Rsd: [20,30) mm','Units','Normalized','Position',...
    [-0.04,(dz3+SepSize+0.6*dz2)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');
h3=text('String','Rsd: [30,40) mm','Units','Normalized','Position',...
    [-0.04,(0.60*dz3)/dzT],'Rotation',90,'Color','w',...
    'FontSize',12,'HorizontalAlignment','center');
end

set(gca,'XTick',[],'YTick',[],'Box','on','Color','w');
colormap(gray(1000))
colorbar('Color','w');
end

set(gcf,'Color',BkndColor)
