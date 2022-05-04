function [Nd, Ns, Nwl, Nt] = Read_AcqDecode_Header(fid)

% READ_ACQDECODE_HEADER Reads the header section of a binary NIRS file.
%
%   [Nd, Ns, Nwl, Nt] = READ_ACQDECODE_HEADER(fid) reads the first bytes of
%   an open binary NIRS file identified by "fid", and selects the correct
%   format to load the rest of the header.
%
% See Also: LOADMULTI_ACQDECODE_DATA, LOAD_ACQDECODE_DATA, FOPEN, FREAD.
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

%% Read header and select correct format.
switch fread(fid, 1, 'uint16') % Read the format string.
    case 44801 % Original
        % 0xAF01
        fseek(fid, 1, 'cof');
        Nd = fread(fid, 1, 'uint16', 2);
        Ns = fread(fid, 1, 'uint16');
        Nwl = 2;
        Nt = [];
    case 44802 % Old
        % 0xAF02
        fseek(fid, 1, 'cof');
        null_bytes = fread(fid, 1, 'uint16', 2);
        Nd = fread(fid, 1, 'uint16', 2);
        Ns = fread(fid, 1, 'uint16');
        Nwl = 2;
        Nt = [];
    case {44803, 44804} % Current
        % 0xAF04
        null_bytes = fread(fid, 1, 'uint16');
        Nd = fread(fid, 1, 'uint16');
        Ns = fread(fid, 1, 'uint16');
        Nwl = fread(fid, 1, 'uint16', 2); % OLD NOTE: the extra bit here has a meaning, but I don't remember it
        Nt = fread(fid, 1, 'uint32');
end



%
