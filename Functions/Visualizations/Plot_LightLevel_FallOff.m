function Plot_LightLevel_FallOff(data,info,params)
%
% This function generates a quick light level fall off. Source-detector
% distances are assumed to be cotained within info.pairs.r2d or
% info.pairs.r3d.

%% Parameters and Initialization
if ~exist('params', 'var')
    params = [];
end

BkgdColor = 'k';
params.fig_size = [200, 200, 560, 420];
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end

Nm = sum(info.pairs.WL==1);
wls=unique(info.pairs.lambda);
Nwls=length(wls);
leg=cell(Nwls,1);
for j=1:Nwls, leg{j}=[num2str(wls(j)),' nm'];end
Phi_0=mean(abs(data),2);
M=ceil(max(log10(abs(data(:)))));
yM=10^M;
if isfield(info.pairs,'r3d')
    if any(info.pairs.r3d)
        r=info.pairs.r3d(info.pairs.lambda==wls(1));
    end
else
    r=info.pairs.r2d(info.pairs.lambda==wls(1));
end

xM=min([max(r)+5,100]);

%% Draw Plot
semilogy(r,reshape(Phi_0,Nm,[]),'.');
axis([0,xM,1e-6,yM])
xlabel('Source-Detector Separation ( mm )','Color','w')
ylabel('{\Phi_0} ( {\mu}W )','Color','w')
set(gca,'XColor','w','YColor','w','Xgrid','on','Ygrid','on','Color','k')
legend(leg,'Color','w')