function [HbO,HbR] = MBLL(data_in, info, ppf, e_mat)
% Convert OD to concentrations using the Modified Beer-Lambert Law.
%
% Takes the optical density data (data_in) measured at two wavelengths for 
% each source-detector pair and multiplies it by the inverse of the 
% extinction coefficient matrix (e_mat). Then, divides by the partial
% path length factor (ppf) multiplied by the source-detector separation.
% Creates the outputs HbO and HbR, the concentrations of 
% oxygenated and deoxygenated hemoglobin, respectively.
%
% INPUTS:
% data_in: optical density data in NeuroDOT format (# channels x time) 
% info: NeuroDOT info structure
% ppf: Partial path length factors for each wavelength. [1x # wavelengths]
%  e_mat: spectroscopy matrix
%
% OUTPUTS:
% HbO: array with concentration of oxy-Hb for each channel
% HbR: array with concentration of deoxy-Hb for each channel
%
%
%% Parameters and initialization
data_in = data_in';
nTpts = size(data_in,1);
Lambda = unique(info.pairs.lambda); 
nWav   = length(Lambda);
lst = find(info.pairs.WL==1);  % Find indices of first wavelength
Nm = length(lst);
HbO = zeros(Nm, nTpts); % Initialize empty HbO and HbR arrays
HbR = HbO;

if length(ppf) < nWav
    warning(['Length of ppf does not match the number of wavelengths.', ...
        'Falling back to ppf=1 for all wavelengths.']);
    ppf = ones(1, nWav);
elseif length(ppf) > nWav
    d = length(ppf)-nWav;
    ppf = ppf(1:d);
end


%% MBLL Algorithm
einv = inv(e_mat'*e_mat)*e_mat'; % Invert extinction coefficient matrix

for idx=1:Nm
    idx1 = lst(idx); % Index of WL1 channel
    idx2 = find(info.pairs.WL>1 & info.pairs.Src==info.pairs.Src(idx1)...
        & info.pairs.Det==info.pairs.Det(idx1)); % Index of corresponding WL2 channel
    if ppf(1)~=1
        temp = (einv * (data_in(:,[idx1, idx2'])./(ones(nTpts,1)*info.pairs.r3d(idx)*ppf))')';
    else
        temp = (einv * (data_in(:,[idx1, idx2']))')';
    end
       HbO(idx,:) = temp(:,1);
       HbR(idx,:) = temp(:,2);
end   


end
