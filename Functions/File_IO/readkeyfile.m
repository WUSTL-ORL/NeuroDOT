function [fvalues]=readkeyfile(filename)

% readkeyfile reads in key-files and exports a struct with the included
% information. Key-files are NeuroDOT's main information storage format to
% talk with each other and to input processing commands into the Matlab
% experimental design programs. Basically, each file contains a list of
% variables and the values to assign to them. readkeyfile assigns each
% variable to a field in a struct and gives that field the requested value.
% 
% Syntax: >> [fvalues]=readkeyfile(filename)
% 
% filename is a struct containing the filename (with path and extension).
% fvalues is a struct containing the returned variables. The fields of this
% struct depend on the keys in the file. The interpretation of various keys
% is handled by the sub-program interpstring. This includes a list of known
% strings and their formatting. Unknown strings are given field names that
% are simply the lower-case of their key name, and they are converted to
% doubles unless they contain any letters, in which case, they are retained
% as strings. 
%
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

fvalues=struct; % initialize

if ~isstr(filename); error('** You entered a filename that was not a string **'); end

% if no extension is given, assume *txt
if ~strcmp(filename((end-3)),'.'); filename=strcat(filename,'.txt'); end

fid=fopen(filename);

if fid==-1
    error(['** The key-file you requested (',filename,') does not exist **'])
end

readflag=1;
while readflag
    readstring=fgetl(fid);
    if readstring==-1 % End of File
        readflag=0;
    elseif (strncmp(readstring,'#',1) || strncmp(readstring,'%',1)) % Comment Line
        continue
    elseif strncmp(readstring,'!',1) % Include additional nested key-file
        fvalues2=readkeyfile(readstring(3:end));
        fvalues=mergestruct(fvalues,fvalues2,'last');
    elseif isempty(readstring) % Blank Line
        continue
    else % if line with KEY VALUE
        [key value]=interpstring(readstring); % interpret line
        fvalues.(key)=value; % enter into struct
    end
end

fclose(fid);

end

%% interpstring()
function [key2,value]=interpstring(rs)

w=isstrprop(rs,'wspace'); % find white-space
kvsplit=find(w,1); % KEY/VALUE divide is the first white-space

key=rs(1:(kvsplit-1)); % Extract KEY
value=rs((kvsplit+1):end); % Extract VALUE (as string)

switch key
    % New naming convetion
    case {'ENCODING_NAME','PAD_NAME'}
        key2=lower(key(1:3));
    case 'NCOLOR'
        key2='cnum';
        value=str2double(value);
    case {'NDET','NSRC'}
        key2=[lower(key(2:end)),'num'];
        value=str2double(value);
    case {'P','W'}
        key2=key;
        value=str2double(value);
        
    % Use same name, keep value as string
    case {'NAME','DATE','TIME','TAG','COMMENT','GI','E','CENT','SMETH',...
            'FLAGDIR','DIR','PATH','BASE','SHELL','LOGIC','INVERT',...
            'START_TIME','END_TIME'}
        key2=lower(key);
        
    % Use same name, change value to number
    % Scalars
    case {'UNIX_TIME','RUN','NREG','NTS','NSAMP','NBLANK','NMOTU','NAUX',...
            'DET','LOWPASS1','LOWPASS2','LOWPASS3','HIGHPASS','GSR','BTHRESH',...
            'OMEGA_LP1','OMEGA_LP2','OMEGA_LP3','OMEGA_HP','FREQOUT','STARTPT',...
            'CTHRESH','SEEDSIZE','GLOBREG','WRITECENT','STHRESH','RSTOL',...
            'PWAVE','PRATE','ABBELT','THBELT','PNCYCLE','HANDTHRESH',...
            'GBOX','GSIGMA','SHELLIN','SHELLOUT','ALIAS','NDIO',...
            'STEPL','BUFFL','FRAMEHZ','BUFFPERCENT','START_UNIX_TIME'}
        if strcmp(key,'START_UNIX_TIME'), key='UNIX_TIME';end
        key2=lower(key);
        value=str2double(value);
    % Vectors and Exponentials
    case {'FREQ','DIV','NN1','NN2','LAMBDA1','LAMBDA2','LEDHZ'}
        key2=lower(key);
        value=str2num(value);
        
    otherwise
        % disp(['** readkeyfile() encountered an unknown KEY: ',key,' **'])      
        key2=lower(key);
        if any(isstrprop(value,'alpha')) % Keep values with letters as strings
            % disp('  Maintaining ASCII format')
        else % Others turn into numbers
            % disp('  Converting to Double format')
            value2=str2num(value);
            if isempty(value2)
                eval(['value2=',value])
            end
            value=value2;
        end
end
end


function [Sout]=mergestruct(varargin)

% mergestruct is a program to combine multiple structs into one. Syntax:
% 
% >> [Sout]=mergestruct(S1,S2,S3,...,option)
% 
% S1, S2, and S3 are input structs. Sout is the output struct. Every field
% from the input structs is included in Sout with the same values as were
% found in the input structs. You can include as many input structs as
% desired.
% 
% mergestruct checks to make sure that no fields are duplicated between the
% inputs. If the same field is present in both input structs, but with
% different values, then option (which can be included anywhere, not just
% at the end) tells mergestruct how to deal with the overlap. Valid values
% are: 'first' (use the value from the first struct in the list), 'last'
% (use the value from the last struct in the list), or 'exclude' (don't
% include this field in the final struct). 
%

% Process variable inputs
numstruct=0; % initialize
for n=1:nargin % loop over all inputs
    switch class(varargin{n})
        case 'struct' % if struct, then put into array of structs.
            numstruct=numstruct+1; % increment counter
            Sin{numstruct}=varargin{n}; % Sin{} will hold all structs
        case 'char' % if string, then interpret command.
            oflag=varargin{n}; % only command is "overlap flag"
    end
end

% default
if ~exist('oflag','var'); oflag='first'; end 

switch numstruct
    case 0 % if not structs were input
        Sout=[];
        return
    case 1 % if only one input struct
        Sout=Sin{1};
        return
    otherwise % actual merging!
        Sout=Sin{1}; % start with the first struct
        lostnum=0; % initialize number of overlaps
        lostlist=[]; % initialize list of overlaps
        for n=2:numstruct % loop over rest of structs
            [Sout lostnum lostlist]=combstruct(Sout,Sin{n},oflag,lostnum,lostlist);
        end
end

% report on overlaps
if lostnum
    disp('** Multiple input structs contain the same fields with different values **')
    disp(lostlist)
end

end

%% combstruct()
% Merge two structs: S2 into Sout
function [Sout lostnum lostlist]=combstruct(Sout,S2,oflag,lostnum,lostlist)

% get new field names
N2=fieldnames(S2);

for n=1:numel(N2) % loop over fields
    fname=N2{n}; % current field
    
    % is this field in Sout already?
    if isfield(Sout,fname) % if it is, we need to check if the values are equal
        [cflag]=checkoverlap(Sout.(fname),S2.(fname)); % check equality
        
        % if no conflict, then we can ignore it
        % value is already in Sout
        
        % if there is a conflict
        if cflag
            lostnum=lostnum+1; % incremement loss number
            lostlist{lostnum}=fname; % add to list
            
            % resolve
            if strcmp(oflag,'exclude') % if we want to exclude
                Sout=rmfield(Sout,fname); % then delete it from the structs
            else % otherwise
                Sout.(fname)=resolveoverlap(Sout.(fname),S2.(fname),oflag);
            end
        end
    else % if no overlap
        Sout.(fname)=S2.(fname); % then we're fine, just add it to Sout
    end
end

end

%% checkoverlap()
% check to see if values V1 and V2 are the same
function [cflag]=checkoverlap(V1,V2)

cflag=0; % initialize
if ~strcmp(class(V1),class(V2)) % if values are of different types (e.g., number vs. string)
    cflag=1; % then, they definitely aren't equal
    return
else % then compare values
    if ischar(V1) % check equality of strings
        if ~strcmp(V1,V2) % compare strings
            cflag=1;
            return
        end
    elseif isnumeric(V1) % check equality of numbers
        if V1~=V2 % compare numbers
            cflag=1;
            return
        end
    end
end

end

%% resolveoverlap()
function Vout=resolveoverlap(V1,V2,oflag)

switch oflag
    case 'first'
        Vout=V1;
    case 'last'
        Vout=V2;
end

end
