function Save_HOMER(data, info, filename, pn)

% SAVE_HOMER Converts an ND2 workspace to HOMER2 ".nirs" format to save.
% 
%   SAVE_HOMER(data, info, filename, pn) converts the ND2 variables "data"
%   and "info" into HOMER2's ".nirs" format, and saves into a file
%   specified by "filename" and path "pn".
% 
% Dependencies: CONVERTER_ND2_TO_HOMER.
% 
% See Also: LOAD_HOMER.
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
if ~exist('pn', 'var')  ||  isempty(pn)
    [pn, filename, ext] = fileparts(filename);
else
    [pn, filename, ext] = fileparts(fullfile(filename, pn));
end

%% Convert to HOMER format.
nirs = converter_ND2_to_HOMER(data, info);

%% Save HOMER .nirs file.
save(fullfile(pn, [filename, ext]), '-struct', 'nirs', '-mat');



%
