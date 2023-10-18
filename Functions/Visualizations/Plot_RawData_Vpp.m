function Plot_RawData_Vpp(data,info,params)
%
% This function generates a figure of the 'Vpp', which are essentially the
% mean light level for each source (1st fig) and each detector
% (2nd fig).
%   data        assumed to be raw data (not optical density) 
%                   in measurements by time
%   info        standard info structure. assumes to contain info.pairs
%                   which contains the relevant measurement list and
%                   info.MEAS which contains a logical array Clipped which
%                   the same size as the measure. If info.MEAS.Clipped is
%                   not present, it is assumed there are no clipped
%                   channels.
%   params      optional input to set flags for visualizations. fields of
%                   note listed here:
%
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

Ns=max(info.pairs.Src);
Nd=max(info.pairs.Det);
Sxm=0;
SxM=Ns+1;
Dxm=0;
DxM=Nd+1;
wls=unique(info.pairs.WL);
Nwls=length(wls);
Cmarks={'b','r','m','c'};
Smarks={'s','o','+','d'};

if isfield(info,'MEAS')
    if isfield(info.MEAS,'Clipped')
        clipped=info.MEAS.Clipped;
    else
        clipped=zeros(Ns*Nd*Nwls,1)==1;
    end
else
    clipped=zeros(Ns*Nd*Nwls,1)==1;
end

if ~exist('params','var'),params=struct;end
if ~isfield(params,'ym'), params.ym=1e-6;end
if ~isfield(params,'yM'), params.yM=1e1;end
if ~isfield(params,'NN'), params.NN=1;end

if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end


%% Calc mean log mean light level
Phi_0=mean(data,2);


%% Plot data
for j=1:Nwls
subplot(2,1,1); % Sources Vpp
keep1=info.pairs.WL==wls(j) & info.pairs.NN==params.NN & ~clipped;
semilogy(info.pairs.Src(keep1),(Phi_0(keep1)),[Cmarks{j},Smarks{j}]);hold on

keep1=info.pairs.WL==wls(j) & info.pairs.NN==params.NN & clipped;
semilogy(info.pairs.Src(keep1),(Phi_0(keep1)),['kx']);

xlabel('Sources');ylabel(['\Phi_0'])
axis([Sxm,SxM,params.ym,params.yM])


subplot(2,1,2); % Detectors Vpp
keep1=info.pairs.WL==wls(j) & info.pairs.NN==params.NN & ~clipped;
semilogy(info.pairs.Det(keep1),(Phi_0(keep1)),[Cmarks{j},Smarks{j}]);hold on

keep1=info.pairs.WL==wls(j) & info.pairs.NN==params.NN & clipped;
semilogy(info.pairs.Det(keep1),(Phi_0(keep1)),['kx']);

xlabel('Detectors');ylabel(['\Phi_0'])
axis([Dxm,DxM,params.ym,params.yM])
end