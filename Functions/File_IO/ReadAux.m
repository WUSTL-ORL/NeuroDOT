function aux = ReadAux(info, filename, pn)

% READAUX Reads auxiliary channel files.
% 
%	aux = READAUX(info, filename, pn) reads auxiliary channel files contained
%	in the -info.txt key specified by "info.io.run", "filename", and path
%	"pn". The channel files must be located in the same directory as the
%	raw data.
% 
% See Also: LOADMULTI_ACQDECODE_DATA, LOAD_ACQDECODE_DATA, READINFOTXT.
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

%% Find all files that match in this folder.
pn_new = pn;

temp = what(pn_new); % Strips out the '\'.

name = fullfile(temp.path,...
    [filename(1:6), '-run', num2str(info.io.run, '%03g'), '*.raw']);

aux_chan = dir(name);

if isempty(aux_chan)
    error(['** ReadAux: No .raw files were found for ', name, ' ***'])
end

if size(aux_chan,1)>1
    aux_chan=aux_chan(size(aux_chan,1)-1);
end

%% Load files.
n = 0;
for chan = aux_chan
    n = n + 1;
    fid = fopen(fullfile(temp.path, chan.name));
    aux.(['aux', num2str(n)]) = fread(fid, 'float32');
    fclose(fid);
end



%
