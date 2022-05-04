function data_out = converter_data(data_in, conversion, info)

% CONVERTER_DATA Reshapes data between ND1 and ND2 formats.
% 
%   data_out = CONVERTER_DATA(data_in, conversion, info) reshapes any ND1 array
%   "data_in" to ND2 array format in "data_out", or vice versa. The
%   direction is determined by the string "conversion", and in the 'ND2 to
%   ND1' direction, "info" is needed.
% 
%   "conversion" can be either 'ND2 to ND1' or 'ND1 to ND2'.
% 
% See Also: CONVERTER_ND1_TO_ND2, CONVERTER_ND2_TO_ND1, CONVERTER_INFO.
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

%% Reshape to new format.
switch conversion
    case 'ND1 to ND2'
        [Nm, Nwl, Nt] = size(data_in);
        if isstruct(data_in)
            data_out.sin = reshape(data_in.sin, [], Nt);
            data_out.cos = reshape(data_in.cos, [], Nt);
        else
            data_out = reshape(data_in, [], Nt);
        end
    case 'ND2 to ND1'
        Nwl = length(unique(info.pairs.WL));
        [Nm, Nt] = size(data_in);
        if isstruct(data_in)
            data_out.sin = reshape(data_in.sin, Nm/2, Nwl, Nt);
            data_out.cos = reshape(data_in.cos, Nm/2, Nwl, Nt);
        else
            data_out = reshape(data_in, Nm/2, Nwl, Nt);
        end
end



%
