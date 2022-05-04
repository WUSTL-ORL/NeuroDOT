function hem = gethem(data, info, sel_type, value)

% GETHEM Calculates mean across set of measurements.
%
%   hem = GETHEM(data, info) takes a light-level array "data" of the format
%   MEAS x TIME, and using the scan metadata in "info.pairs" averages the
%   measurements for each wavelength present. The result is
%   referred to as the "hem" of a measurement set. If there is a
%   good measurements logical vector present in "info.MEAS.GI", it will be
%   applied to the data; otherwise, "info.MEAS.GI" will be set to true for
%   all measurements (i.e., all measurements are assumed to be good). The
%   variable "hem" is output in the format WL x TIME.
% 
%   hem = GETHEM(data, info, sel_type, value) allows the user to set the
%   criteria for determining shallow measurements. "sel_type" can be 'r2d',
%   'r3d', or 'NN', corresponding to the columns of the "info.pairs" table,
%   and "value" can either take the form of a two-element "[min, max]"
%   vector (for 'r2d' and 'r3d'), or a scalar or vector containing all
%   nearest neighbor numbers to be averaged. By default, this function
%   averages the first nearest neighbor.
%
% See Also: REGCORR, DETREND_TTS.
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

%% Parameters and Initialization.
[Nm, Nt] = size(data);
cs = unique(info.pairs.WL); % WLs.
Nc = length(cs); % Number of WLs.
hem = zeros(Nc, Nt);

if exist('sel_type', 'var')  &&  exist('value', 'var')
    if ischar(value)
        value = str2double(value);
    end
else
    sel_type = 'r2d';
    value = [10, 16];
end

switch sel_type
    case 'r2d'
        keep_R_NN = (info.pairs.r2d >= value(1))  &...
            (info.pairs.r2d <= value(2));
    case 'r3d'
        keep_R_NN = (info.pairs.r3d >= value(1))  &...
            (info.pairs.r3d <= value(2));
    case 'NN'
        keep_R_NN = zeros(size(info.pairs.NN));
        for k = value
            keep_R_NN = keep_R_NN  |  (info.pairs.NN == k);
        end
end

%% If no GI, assume all measurements good.
if (isfield(info, 'MEAS')  &&  ~isfield(info.MEAS, 'GI'))
    info.MEAS.GI = true(Nm, 1);
elseif ~isfield(info, 'MEAS')
    info.MEAS = struct;
    info.MEAS.GI = true(Nm, 1);
end

%% Make "hem" by averaging measurments (e.g., nn1).
for k = 1:Nc
    keep = keep_R_NN  &  info.pairs.WL == cs(k)  &  info.MEAS.GI;
    hem(k, :) = mean(data(keep, :), 1);
end



%
