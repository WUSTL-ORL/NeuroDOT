function io = ReadInfoTxt(filename)

% READINFOTXT Reads a key file.
% 
%   io = READINFOTXT(filename) scans the -info.txt key file specified by
%   "filename" and parses its contents into key-value pairs, which are
%   returned in the "io" substructure of the main ND2 "info" structure.
%   Most of the pairs are either only used in the loading process or kept
%   for later troubleshooting if the data set shows inconsistencies.
% 
% See Also: LOADMULTI_ACQDECODE_DATA, LOAD_ACQDECODE_DATA, READAUX.
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

%% Parameters and initialization.
io = struct;

%% Read text file.
fid = fopen(filename);
temp = textscan(fid, '%s %s');
fclose(fid);
numKey=0;

%% Pairs each key-value pair.
for pairs = [temp{1}'; temp{2}']
    % Note: I would add a third row of proper info fields, but the order of
    % the pairs is indeterminate, so we couldn't sync such a row to the
    % rest of the pairs.
    
    key = pairs{1};
    val = pairs{2};
    
    switch key
        case {'ENCODING_NAME', 'PAD_NAME'}
            io.(lower(key(1:3))) = val;
        case 'NCOLOR'
            io.Nwl = str2double(val);
        case 'NDET'
            io.Nd = str2double(val);
        case 'NSRC'
            io.Ns = str2double(val);
        case 'ENCPASSN'
            io.EncPassN = str2double(val);
        case {'P', 'W'}
            io.(key) = str2double(val);
        case {'NAME', 'DATE', 'TIME', 'TAG', 'COMMENT', 'GI', 'E', 'CENT', 'SMETH',...
                'FLAGDIR', 'DIR', 'PATH', 'BASE', 'SHELL', 'LOGIC', 'INVERT',...
                'START_TIME', 'END_TIME'}
            % Use same name, keep value as string
            io.(lower(key)) = val;
        case {'UNIX_TIME', 'RUN', 'NREG', 'NTS', 'NSAMP', 'NBLANK', 'NMOTU', 'NAUX',...
                'DET', 'LOWPASS1', 'LOWPASS2', 'LOWPASS3', 'HIGHPASS', 'GSR', 'BTHRESH',...
                'OMEGA_LP1', 'OMEGA_LP2', 'OMEGA_LP3', 'OMEGA_HP', 'FREQOUT', 'STARTPT',...
                'CTHRESH', 'SEEDSIZE', 'GLOBREG', 'WRITECENT', 'STHRESH', 'RSTOL',...
                'PWAVE', 'PRATE', 'ABBELT', 'THBELT', 'PNCYCLE', 'HANDTHRESH',...
                'GBOX', 'GSIGMA', 'SHELLIN', 'SHELLOUT', 'ALIAS', 'NDIO',...
                'STEPL', 'BUFFL', 'FRAMEHZ', 'BUFFPERCENT', 'START_UNIX_TIME'}
            % Scalars, use same name, change value to number
            if strcmp(key, 'START_UNIX_TIME')
                key = 'UNIX_TIME';
            end
            io.(lower(key)) = str2double(val);
        case {'FREQ', 'DIV', 'NN1', 'NN2', 'LAMBDA1', 'LAMBDA2', 'LEDHZ'}
            % Vectors and Exponentials, same as above.
            io.(lower(key)) = str2double(val);
        case {'MARK'} % --> grab all remaining fields and put in as string
            MarkIdx=find(ismember(temp{1,1},'MARK'));
            msg=[];
            maxIdx=size(temp{1,1},1);
            for j=MarkIdx:maxIdx
                msg=strcat(msg,temp{1,1}(j),'_',temp{1,2}(j),'_');
            end
            io.(lower(key)) =msg;
            return
        otherwise
            if ~any(isstrprop(val, 'alpha')) % Keep values with letters as strings
                val2 = str2double(val);
                if isempty(val2)
                    eval(['value2 = ', val]);
                end
                val = val2;
            end
            io.(lower(key)) = val;
    end
end



%
