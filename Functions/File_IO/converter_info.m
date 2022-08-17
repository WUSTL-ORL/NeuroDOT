function info_out = converter_info(info_in, conversion)

% CONVERTER_INFO Converts between ND1 and ND2 info structure formats.
%
%   info_out = CONVERTER_INFO(info_in, conversion) adapts the metadata
%   structures between ND1 and ND2.
%
%   "conversion" can be either 'ND2 to ND1' or 'ND1 to ND2'.
% 
%   All leftover fields in 'ND1 to ND2' are stored in "info_out.misc".
%
%   Only intended for use with ND1 and ND2 data structures.
%
% See Also: CONVERTER_DATA, CONVERTER_ND1_TO_ND2, CONVERTER_ND2_TO_ND1.
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

%% Get all fieldnames in "info".
inames = fieldnames(info_in);
wl1_default = 750;
wl2_default = 850;

info_out = [];
switch conversion
    case 'ND1 to ND2'
        %% Build "info.paradigm".
        WANTparadigm = {'synchpts', 'synchtype', 'Pulse_1', 'Pulse_2',...
            'Pulse_3', 'Pulse_4'};
        NEWparadigm = WANTparadigm;
        for str = [WANTparadigm; NEWparadigm]
            tf = strcmp(str{1}, inames);
            if any(tf)
                info_out.paradigm.(str{2}) = info_in.(str{1});
                inames(tf) = [];
            end
        end
        if isfield(info_out, 'paradigm')  &&  isfield(info_out.paradigm,...
                'synchpts')
            info_out.paradigm.init_synchpts = info_out.paradigm.synchpts;
        end
        
        %% Build "info.system".
        WANTsys = {'pad', 'framerate', 'currHz'}; %#ok<*CCAT>
        NEWsys = {'PadName', 'init_framerate', 'framerate'};
        for str = [WANTsys; NEWsys]
            tf = strcmp(str{1}, inames);
            if any(tf)
                info_out.system.(str{2}) = info_in.(str{1});
                inames(tf) = [];
            end
        end
        if ~isfield(info_out.system,'framerate')
            info_out.system.framerate=info_out.system.init_framerate;
        end
        
    if ~isfield(info_in,'Pad_info')
        %% Build "info.optodes".
        if isfield(info_in, 'pad')
            info_out.optodes.CapName = info_in.pad;
        end
        WANToptodes = {'dpos', 'dpos3', 'spos', 'spos3'};
        NEWoptodes = {'dpos2', 'dpos3', 'spos2', 'spos3'};
        if isfield(info_in, 'grid')  &&  ~isempty(info_in.grid)
            gnames = fieldnames(info_in.grid);
            for str = [WANToptodes; NEWoptodes]
                if any(strcmp(str{1}, gnames))
                    info_out.optodes.(str{2}) = info_in.grid.(str{1});
                end
            end
            tf = strcmp('grid', inames);
            inames(tf) = [];
        end
        info_out.optodes.plot3orientation.i = 'R2L';
        info_out.optodes.plot3orientation.j = 'P2A';
        info_out.optodes.plot3orientation.k = 'D2V';
        
        %% Build "info.pairs" table.
        info_out.pairs = struct;
        
        % Check for iRad and its subfields. Cannot do this without that!
        rname = 'iRad';
        tf = strcmp(rname, inames);
        if ~any(tf)
            tf = strcmp('Rad', inames);
            rname = 'Rad';
        end
        if any(tf)
            WANTirad = {'meas', 'nn1', 'nn2', 'nn3', 'nn4', 'nn5'};
            irnames = fieldnames(info_in.(rname));
            tf2 = [];
            for str = WANTirad
                tf2(end+1, :) = strcmp(str{1}, irnames);
            end
            if all(any(tf2'))
                
                info_out.pairs.Src = info_in.(rname).meas(:, 1);
                info_out.pairs.Det = info_in.(rname).meas(:, 2);
                
                Nm = length(info_in.(rname).meas);
                
                info_out.pairs.NN = zeros(Nm, 1);
                for k = 1:5
                    info_out.pairs.NN(info_in.(rname).(['nn' num2str(k)])) = k;
                end
                info_out.pairs = [info_out.pairs; info_out.pairs];
                
                % This section assumes WLs of 750 and 850, and a continuous 
                % wave modulation scheme. Can be changed if necessary.
                info_out.pairs.WL(1:Nm, 1) = 1;
                info_out.pairs.WL(Nm+1:Nm*2, 1) = 2;
                info_out.pairs.lambda(1:Nm, 1) = wl1_default;
                info_out.pairs.lambda(Nm+1:Nm*2, 1) = wl2_default;
                info_out.pairs.Mod(1:Nm*2, 1) = {'CW'};
                
                for k = 1:Nm*2
                    for d = 2:3
                        dd = num2str(d);
                        if isfield(info_out.optodes, ['spos', dd])...
                                &&  isfield(info_out.optodes, ['dpos', dd])
                            info_out.pairs.(['r', dd, 'd'])(k) =...
                                norm([info_out.optodes.(['spos', dd])(info_out.pairs.Src(k),:)-...
                                 info_out.optodes.(['dpos', dd])(info_out.pairs.Det(k),:)]);
                        end
                    end
                end
            end
            info_out.misc.iRad = info_in.(rname);
        end
        
    else
        info_out.pairs=info_in.Pad_info.pairs;
        info_out.optodes=info_in.Pad_info.optodes;
        
    end
        
        %% Build "info.MEAS" if good indices present.
        tf = strcmp('goodindex', inames);
        if any(tf)
            H = height(info_out.pairs);
            info_out.MEAS = struct;
            info_out.MEAS.GI(info_in.goodindex.c1) = true;
            info_out.MEAS.GI(info_in.goodindex.c2 + H/2) = true;
            inames(tf) = [];
        end
        
        %% Grab "info.t4" or create it if not included.
        tf = strcmp('t4', inames);
        if any(tf)
            info_out.tissue.affine = info_in.t4;
            inames(tf) = [];
        else
            info_out.tissue.affine = eye(4);
        end
        
        %% Populate info.io.
        a2z = cellstr(char(97:99)')'; % This gets us 'a' through 'c'. Can be expanded later.
        WANTio = a2z;
        NEWio = a2z;
        for str = a2z
            tf = strcmp(str{1}, inames);
            if any(tf)
                info_out.io.(str{1}) = info_in.(str{1});
                inames(tf) = [];
            end
        end
        
        %% Populate "info.misc" with everything not used thus far.
        if ~isempty(inames)
            for str = inames'
                info_out.misc.(str{1}) = info_in.(str{1});
            end
        end
        
        %% Order all fields.
        info_out = orderfields(info_out);
        if isfield(info_out, 'paradigm')
            info_out.paradigm = orderfields(info_out.paradigm);
        end
        if isfield(info_out, 'system')
            info_out.system = orderfields(info_out.system);
        end
        if isfield(info_out, 'optodes')
            info_out.optodes = orderfields(info_out.optodes);
        end
        
    case 'ND2 to ND1'
        %% Load info.io into info2.a and info2.b.
        WANTio = {'Rad', 'date', 'subj', 'suf'};
        NEWio = {'Rad', 'date', 'subj', 'suf'};
        tf = strcmp('io', inames);
        if any(tf)
            if isfield(info_in.io, 'a')
                info_out.a = info_in.io.a;
                for str = [WANTio; NEWio]
                    if isfield(info_in.io.a, str{1})
                        info_out.(str{2}) = info_in.io.a.(str{1});
                    end
                end
            end
            if isfield(info_in.io, 'b')
                info_out.b = info_in.io.b;
            end
            if ~isfield(info_in.io, {'a', 'b'})
                io_names = fieldnames(info_in.io);
                for str = [WANTio; NEWio]
                    if isfield(info_in.io, str{1})
                        info_out.(str{2}) = info_in.io.(str{1});
                    end
                end
            end
            inames(tf) = [];
        end
        
        %% Load info.misc into info2.
        tf = strcmp('misc', inames);
        if any(tf)
            misc_names = fieldnames(info_in.misc);
            for str = misc_names'
                info_out.(str{1}) = info_in.misc.(str{1});
            end
            inames(tf) = [];
        end
        
        %% Load info.system into info2.
        a2z = cellstr(char(97:99)')'; % This gets us 'a' through 'c'. Can be expanded later.
        WANTsys = {'PadName', 'init_framerate', 'framerate'};
        NEWsys = {'pad', 'framerate', 'currHz'};
        tf = strcmp('system', inames);
        if any(tf)
            for str = [WANTsys; NEWsys]
                if isfield(info_in.system, str{1})
                    info_out.(str{2}) = info_in.system.(str{1});
                end
            end
            inames(tf) = [];
        end
        
        %% Load info.paradigm into info2.
        WANTpar = {'synchpts', 'synchtype', 'Pulse_1', 'Pulse_2',...
            'Pulse_3', 'Pulse_4'};
        NEWpar = {'synchpts', 'synchtype', 'Pulse_1', 'Pulse_2',...
            'Pulse_3', 'Pulse_4'};
        tf = strcmp('paradigm', inames);
        if any(tf)
            par_names = fieldnames(info_in.paradigm);
            for str = [WANTpar; NEWpar];
                if any(strcmp(str{1}, par_names))
                    info_out.(str{2}) = info_in.paradigm.(str{1});
                end
            end
            inames(tf) = [];
        end
        
        %% Create info2.goodindex from info.meas.
        % NOTE: Must be done before iRad, which pares down info.pairs.
        tf = strcmp('MEAS', inames);
        if any(tf)  &&  isfield(info_in.MEAS, 'GI')
            info_out.goodindex.c1 = find((info_in.MEAS.GI == 1)  &  (info_in.pairs.WL == 1));
            info_out.goodindex.c2 = find((info_in.MEAS.GI == 1)  &  (info_in.pairs.WL == 2)) - size(info_in.pairs,1);
            inames(tf) = [];
        end
        
        %% Create info2.iRad from info.pairs.
        tf = strcmp('pairs', inames);
        if any(tf)
            if ~isfield(info_out, 'iRad')
                info_in.pairs(height(info_in.pairs)/2+1:end, :) = [];
                
                for k = 1:5
                    info_out.iRad.(['nn' num2str(k)]) = find(info_in.pairs.NN == k)';
                end
                
                info_out.iRad.nn12 = sort([info_out.iRad.nn1, info_out.iRad.nn2]);
                info_out.iRad.nn123 = sort([info_out.iRad.nn1, info_out.iRad.nn2,...
                    info_out.iRad.nn3]);
                info_out.iRad.nn1234 = sort([info_out.iRad.nn1, info_out.iRad.nn2,...
                    info_out.iRad.nn3, info_out.iRad.nn4]);
                info_out.iRad.nn12345 = sort([info_out.iRad.nn1, info_out.iRad.nn2,...
                    info_out.iRad.nn3, info_out.iRad.nn4, info_out.iRad.nn5]);
                
                info_out.iRad.meas(:, 1) = info_in.pairs.Src;
                info_out.iRad.meas(:, 2) = info_in.pairs.Det;
                
                if isfield(info_in.pairs, 'r2d')
                    info_out.iRad.ir = info_in.pairs.r2d;
                    info_out.iRad.r = info_out.iRad.ir;
                end
            end
            
            inames(tf) = [];
        end
        
        %% Create info2.grid from info.optodes.
        WANToptodes = {'dpos2', 'dpos3', 'spos2', 'spos3'};
        NEWoptodes = {'dpos', 'dpos3', 'spos', 'spos3'};
        tf = strcmp('optodes', inames);
        if any(tf)
            gnames = fieldnames(info_in.optodes);
            for str = [WANToptodes; NEWoptodes]
                if any(strcmp(str{1}, gnames))
                    info_out.grid.(str{2}) = info_in.optodes.(str{1});
                end
            end
            inames(tf) = [];
        end
        % Get detector number.
        poscel = {'dpos', 'dpos3'};
        tf = isfield(info_out.grid, poscel);
        if any(tf)
            idx = find(tf == 1, 1, 'first');
            info_out.grid.detnum = size(info_out.grid.(poscel{idx}), 1);
        end
        % Get source number.
        poscel = {'spos', 'spos3'};
        tf = isfield(info_out.grid, poscel);
        if any(tf)
            idx = find(tf == 1, 1, 'first');
            info_out.grid.srcnum = size(info_out.grid.(poscel{idx}), 1);
        end
        
        %% Turn info.tissue.affine into info2.t4.
        tf = strcmp('tissue', inames);
        if isfield(info_in, 'tissue')
            if isfield(info_in.tissue, 'affine')
                info_out.t4 = info_in.tissue.affine;
            elseif isfield(info_in.tissue, 't4') % Backw
                info_out.t4
            end
            inames(tf) = [];
        end
        
        %% Anything that's not already in info2 goes in now.
        if ~isempty(inames)
            for str = inames'
                info_out.(str{1}) = info_in.(str{1});
            end
        end
        
        %% Reorder fields.
        info_out = orderfields(info_out);
        
        
end



%
