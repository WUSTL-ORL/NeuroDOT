function [data, info, aux] = converter_HOMER_to_ND2(nirs)

% CONVERTER_HOMER_TO_ND2 Converts workspace from HOMER to ND2.
%
%   [data, info, aux] = CONVERTER_HOMER_TO_ND2(nirs) takes the "nirs"
%   structure loaded from a HOMER ".nirs" file and converts its fields into
%   the ND2 variables "data", "info", and optionally "aux".
% 
% See Also: CONVERTER_ND2_TO_HOMER.
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
data = [];
info = [];

%% Convert "data" first.
try
    assert(isfield(nirs, 'd'))
    data  = nirs.d';
    
    if isfield(nirs, 'tIncMan')
        data = data(:, nirs.tIncMan == 1);
    end
    
    %% Save original structure.
    info.io.nirs = nirs;
    
    %% Get framerate.
    assert(isfield(nirs, 't')) % assert for essential info, isfield for nonessential info.
    info.system.framerate = 1 / mean(diff(nirs.t));
    info.system.init_framerate = info.system.framerate;
    
    %% Copy Good-Meas data.
    if isfield(nirs.SD, 'MeasListAct')
        info.MEAS = table(nirs.SD.MeasListAct, 'VariableNames', {'GI'});
    end
    
    %% Detect which measurement list is present.
    if isfield(nirs, 'MeasList')
        MM = nirs.MeasList;
    elseif isfield(nirs, 'ml')
        MM = nirs.ml;
    else
        assert(false)
    end
    info.pairs = table(MM(:, 1), MM(:, 2), 'VariableNames', {'Src', 'Det'});
    Nm = size(info.pairs, 1);
    
    %% Get optode positions.
    assert(all([isfield(nirs, 'SD'), isfield(nirs.SD, {'SrcPos', 'DetPos'})]))
    info.optodes.spos3 = nirs.SD.SrcPos;
    info.optodes.dpos3 = nirs.SD.DetPos;
    
    %% Read synch array "s".
    assert(isfield(nirs, 's'))
    info.paradigm.synchpts = find(sum(nirs.s, 2) ~= 0);
    if isfield(nirs, 's0')
        info.paradigm.init_synchpts = find(sum(nirs.s0, 2) ~= 0);
    else
        info.paradigm.init_synchpts = info.paradigm.synchpts;
        info.io.nirs.s0 = nirs.s; % Save original "s" since it might be modified by user in ND2.
    end
    info.paradigm.synchtype = zeros(size(info.paradigm.synchpts));
    for k = 1:size(nirs.s, 2)
        % WARNING: If stims are on same time point, this overwrites lower
        % stims. However, we already do not use "synchtype" in most
        % functions, so this should not be an issue most of the time.
        % Instead, USE the "Pulse_[x]" fields, which are complete and allow
        % for overlapping stimulus groups.
        info.paradigm.synchtype = nirs.s(info.paradigm.synchpts, k);
        info.paradigm.(['Pulse_', num2str(k)]) = find(info.paradigm.synchtype == k);
    end
    
    %% Calculate NN's from positions. (note, this will sometimes be inaccurate)
    % First, calculate S-D distances.
    for k = 1:Nm
        r3d(k) = norm(info.optodes.spos3(info.pairs.Src(k), :) - info.optodes.dpos3(info.pairs.Det(k), :));
    end
    radii = unique(r3d);
    % Arbitrary logical test - should probably have no more than 9 nn's.
    if numel(radii) > 9
        % Perform a "round by first" algorithm. IE, (1) take a histogram
        % with 50 bins (so that we achieve a good rounding resolution).
        % (2) Select peaks in histogram. (3) Divide raw data by the first
        % peak, call those NN1. (4) Remove those from the group, then start
        % over with a new histogram of the remaining.
        A = r3d;
        k = 1;
        NN = zeros(Nm, 1);
        while ~all(isnan(A))  &&  (k <= 9)
            [counts, centers] = hist(A, 50);
            [~, locs] = findpeaks(counts, 'MinPeakProminence', 10);
            if isempty(locs)
                NN(~isnan(A)) = k - 1;
                break
            end
            locs = sort(locs);
            nn_round = centers(locs(1));
            rounded_A = round(A / nn_round); % Anything near the first peak should be rounded to 1.
            rounded_A(rounded_A == 0) = 1; % If there are outliers below, we round them up to 1 as well.
            NN(rounded_A == 1) = k;
            A(rounded_A == 1) = NaN;
            A = A - nn_round;
            k = k + 1;
        end
        % This algorithm is sensitive to caps with high pair radius
        % variance between NN groups. This algorithm also is primarily
        % intended for high-density caps. Suggest we recommend to user that
        % they check entries personally, or provide their own cap geometry.
    else
        % If not, assume they're all even.
        NN = zeros(Nm, 1);
        for k = 1:numel(radii)
            NN(r3d == radii(k)) = k;
        end
    end
    
    %% Assemble rest of "info.pairs".
    info.pairs.NN(1:Nm, 1) = NN;
    info.pairs.WL(1:Nm, 1) = MM(:, 4);
    info.pairs.lambda(1:Nm, 1) = nirs.SD.Lambda(MM(:, 4));
    info.pairs.Mod(1:Nm, 1) = {'CW'};
    info.pairs.r2d(1:Nm, 1) = 0;
    info.pairs.r3d(1:Nm, 1) = r3d;
    
    %% Fill in remaining fields with defaults.
    info.optodes.CapName = '';
    info.optodes.spos2 = info.optodes.spos3;
    info.optodes.dpos2 = info.optodes.dpos3;
    info.optodes.plot3orientation.i = 'R2L';
    info.optodes.plot3orientation.j = 'P2A';
    info.optodes.plot3orientation.k = 'D2V';
    
    info.system.PadName = '';
    
    info.tissue.affine = eye(4);
    
    %% Extra HOMER data.
    aux = [];
    if isfield(nirs, 'aux') % Just this (for now). May add more later.
        for k = 1:size(nirs.aux, 2)
            aux.(['aux' num2str(k)]) = nirs.aux(:, k);
        end
    end
    
    %% Order fields.
    info = orderfields(info);
    info.system = orderfields(info.system);
    info.paradigm = orderfields(info.paradigm);
    info.optodes = orderfields(info.optodes);
    
catch err
    rethrow(err)
%     error('*** ''.nirs'' file is missing key variables; loading cannot be completed. ***')
end



%
