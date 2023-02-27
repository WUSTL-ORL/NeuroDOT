function [data, info, synch, aux] = LoadMulti_AcqDecode_Data(filename, pn, flags)

% LOADMULTI_ACQDECODE_DATA Loads and combines data from multiple AcqDecode
% files for a single scan.
%
%   [data, info, synch, aux] = LOADMULTI_ACQDECODE_DATA(filename, pn)
%   finds any acquisitions "filename" located in directory "pn", and loads
%   the raw data files and any other relevant scan files contained therein.
%   The information is returned as a raw light level array "data", a
%   metadata structure "info", stimulus synchronization information
%   "synch", and raw auxiliary file data "aux".
% 
%   It is assumed that all files use the AcqDecode format:
%   "DATE-SUBJECT-TAGletter", with a trailing lowercase letter for multiple
%   acquisitions. Also, separate acquisitions are stored in subfolders
%   named by "DATEletter", EG "150115a" and "150115b".
%
%   Currently supports ".mag" and ".iq" files.
% 
%   [data, info, synch, aux] = LOADMULTI_ACQDECODE_DATA(filename, pn,
%   flags) uses the "flags" structure to specify loading parameters.
% 
%   "flags" fields that apply to this function (and their defaults):
%       Nsys            2           Number of acquisitions.
%
% Dependencies: LOAD_ACQDECODE_DATA, CROP2SYNCH. 
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
max_data_trim = 10; % Corresponds to 1 second.
here = pwd;

if ~exist('pn', 'var')  ||  isempty(pn)
    pn = here;
end

if exist('flags', 'var')
    if ~isfield(flags, 'Nsys')
        flags.Nsys = 2;
    end
else
    flags.Nsys = 2;
end

[~, fn, ext] = fileparts(filename);
if isempty(ext)
    ext = '.mag';
end
if ~isempty(fn)
    cel = strsplit(fn, '-');
    if (numel(cel{1}) == 6)  &&  ~any(isletter(cel{1}))
        scan_date = cel{1};
    end
    
end

%% Load a variable number of data sets.
switch flags.Nsys
    case 1
        [data, info, synch, aux, framepts] =...
            Load_AcqDecode_Data(fullfile(pn, [fn, ext]));
        
        [data, info] = Crop2Synch(data, info, flags);
        
    case 2
        %% Load in data set A
        [data_a, info_a, synch.a, aux.a, framepts.a] =...
            Load_AcqDecode_Data(fullfile(pn, [scan_date, 'a'],...
            [fn, 'a', ext]));
        
        % Crop set to outermost synch pulses.
        [data_a, info_a] = Crop2Synch(data_a, info_a, flags);
        
        %% Load in data set B
        [data_b, info_b, synch.b, aux.b, framepts.b] =...
            Load_AcqDecode_Data(fullfile(pn, [scan_date, 'b'],...
            [fn, 'b', ext]));
        
        % Crop set to outermost synch pulses.
        [data_b, info_b] = Crop2Synch(data_b, info_b, flags);
        
        %% Combine data sets.
        % Adopt the info_a as our baseline.
        info = info_a;
        
        % Store io structures for troubleshooting.
        info.io = [];
        info.io.a = info_a.io;
        info.io.b = info_b.io;
        
        % Trim the sets to the same number of frames (to within a second).
        Nts = [size(data_a, 2), size(data_b, 2)];
        L = min(Nts);
        dNt = max(Nts) - L;
        if dNt <= max_data_trim
            data_a = data_a(:, 1:L);
            data_b = data_b(:, 1:L);
        else
            error(['** The decoded data from the systems is inconsistent by ',...
                num2str(dNt), ' too many frames **'])
        end
        
        % Reshape, merge, check and fix meas-info match, reshape.
        data = [reshape(data_a, [], info_a.io.Nwl, L);...
            reshape(data_b, [], info_b.io.Nwl, L)];
        
        % %% Check data measurement number compared to info file
        if size(data,1)~=size(info.pairs,1)
            InfoList=unique([info.pairs.Src,info.pairs.Det],'rows','stable');            
            Nd=size(info.optodes.dpos3,1);
            Ns=size(info.optodes.spos3,1);
            Dlist=[];
            for j=1:Nd
                Dlist=cat(1,Dlist,[(1:Ns)',repmat(j,[Ns,1])]);
            end
            [Ia,Ib]=ismember(Dlist,InfoList,'rows');Ib(Ib==0)=[];
            data=data(Ia,:,:);
        end
        data = reshape(data, [], L);
        
    case 3
        %% Load in data set A
        [data_a, info_a, synch.a, aux.a, framepts.a] =...
            Load_AcqDecode_Data(fullfile(pn, [scan_date, 'a'],...
            [fn, 'a', ext]));
        
        % Crop set to outermost synch pulses.
        [data_a, info_a] = Crop2Synch(data_a, info_a, flags);
        
        %% Load in data set B
        [data_b, info_b, synch.b, aux.b, framepts.b] =...
            Load_AcqDecode_Data(fullfile(pn, [scan_date, 'b'],...
            [fn, 'b', ext]));
        
        % Crop set to outermost synch pulses.
        [data_b, info_b] = Crop2Synch(data_b, info_b, flags);
        
        %% Load in data set C
        [data_c, info_c, synch.c, aux.c, framepts.c] =...
            Load_AcqDecode_Data(fullfile(pn, [scan_date, 'c'],...
            [fn, 'c', ext]));
        
        % Crop set to outermost synch pulses.
        [data_c, info_c] = Crop2Synch(data_c, info_c, flags);
        
        %% Combine data sets.
        % Adopt the info_a as our baseline.
        info = info_a;
        
        % Store io structures for troubleshooting.
        info.io = [];
        info.io.a = info_a.io;
        info.io.b = info_b.io;
        info.io.c = info_c.io;
        
        % Trim the sets to the same number of frames (to within a second).
        Nts = [size(data_a, 2), size(data_b, 2), size(data_c, 2)];
        L = min(Nts);
        dNt = max(Nts) - L;
        if dNt <= max_data_trim
            data_a = data_a(:, 1:L);
            data_b = data_b(:, 1:L);
            data_c = data_c(:, 1:L);
        else
            error(['** The decoded data from the systems is inconsistent by ',...
                num2str(dNt), ' too many frames **'])
        end
        
        % Reshape, merge, reshape - this matches the data to ND2's
        % "info.pairs" measurement list.
        data = [reshape(data_a, [], info_a.io.Nwl, L);...
            reshape(data_b, [], info_b.io.Nwl, L);...
            reshape(data_c, [], info_c.io.Nwl, L)];
        data = reshape(data, [], L);
        
        
end
if istable(info.pairs)
    info.pairs = table2struct(info.pairs, 'ToScalar', true);
end



%
