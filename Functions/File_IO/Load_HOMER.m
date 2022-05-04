function [data, info] = Load_HOMER(filename, pn)

% LOAD_HOMER Converts and loads a HOMER2 ".nirs" file into ND2.
% 
%   [data, info] = LOAD_HOMER(filename, pn) loads a HOMER2 file specified
%   by "filename" and path "pn", and converts it into the ND2 format output
%   as "data" and "info".
% 
% Dependencies: CONVERTER_HOMER_TO_ND2.
% 
% See Also: SAVE_HOMER.
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

if ~exist('pn', 'var')  ||  isempty(pn)
    [pn, filename, ext] = fileparts(filename);
else
    [pn, filename, ext] = fileparts(fullfile(filename, pn));
end

%% Load HOMER .nirs file.
nirs = load(fullfile(pn, [filename, ext]), '-mat');
% We want to load it into a structure, not into the workspace directly.

%% Convert to ND2 format.
[data, info] = converter_HOMER_to_ND2(nirs);



%
