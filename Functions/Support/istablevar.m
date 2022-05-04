function tf = istablevar(table, varname)

% ISTABLEVAR Checks if a variable or list of variables exist in a table.
% 
%   ISTABLEVAR(table, varname) performs the same logical check on tables as
%   ISFIELD does for structure arrays, by checking if the variable(s)
%   "varname" exist within the table "table". "varname" must be either a
%   string or a 1D cell array of strings. A horizontal cell array will give
%   a logical array for each cell array element, while a vertical array
%   will be equivalent to "ALL(ISTABLEVAR)".
% 
%   Example: If info.iRad is a table with variable names 'Src', 'Det',
%   'WL', 'NN', and 'r', then,
% 
%       ISTABLEVAR(info.iRad, {'WL'})
%       
%       gives
%       
%       ans = 
%           1
% 
%       and with a horizontal cell array,
% 
%       ISTABLEVAR(info.iRad, {'WL', 'NN'})
% 
%       gives
% 
%       ans = 
%           1   1
% 
%       However, note that the component function STRCMP in ISTABLEVAR is
%       case-sensitive, so,
% 
%       ISTABLEVAR(info.iRad, {'wL', 'Nn'})
% 
%       will give 
%
%       ans =
%           0   0
% 
%       and with a vertical input array will result in a single output,
%       only returning "true" if both strings are matches.
% 
%       ISTABLEVAR(info.iRad, {'WL'; 'NN'})
% 
%       ans =
%           1
%   Because structures are indexed like tables, if the input is a structure 
%   and not a table, the function will return true if the variable names
%   are fields in the structure of table.
%
%
% See Also: IS*, ISFIELD.
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

if ~istable(table)
    if isstruct (table)
        tf=isfield(table,varname);
        return
    end
end

if ~ischar(varname)  &&  ~iscellstr(varname)
    error('Input "varname" is not a cell array of strings.')
end

if ischar(varname)
    tf = istablevar(table, {varname});
    return
end

tf = [];
for var = varname
    tf(end+1) = any(strcmp(varname{1}, table.Properties.VariableNames));
end



%
