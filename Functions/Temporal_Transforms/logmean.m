function [data_out,Phi_0] = logmean(data_in)

% LOGMEAN Takes the log-ratio of raw intensity data.
%
%   data_out = LOGMEAN(data_in) takes a light-level data array "data_in" of
%   the format MEAS x TIME, and takes the negative log of each element of a
%   row divided by that row's average. The result is output into "data_out"
%   in the same MEAS x TIME format.
%
%   The formal equation for the LOGMEAN operation is:
%        Y_{out} = -log(phi_{in} / <phi_{in}>)
%
%   If the raw optical data phi is complex (as in the frequency domain 
%   case), Y behaves a bit differently. Phi can be defined in terms of Real
%   and Imaginary parts: Phi = Re(Phi) + 1i.*Im(Phi), or in terms of it's
%   magnitude (A) and phase (theta): Phi = A.*exp(1i.*theta).
%   The temporal average of Phi (what we use for baseline) is best
%   calculated on the Real/Imaginary decription: 
%
%       Phi_0=<phi>=mean(data_in,2) = A_o*exp(i*(th_o));
%
%   Taking the logarithm of complex ratio:
%
%       Y_Rytov=-log(phi/<phi>)=-log[A*exp(i*th)/A_o*exp(i*th_o)]
%                              =-[log(A/A_o) + i(th-th_o)];
%
%       Y_Rytov_Re=-log(abs(data_in/Phi_0));
%       Y_Rytov_Im=-angle(data_in/Phi_0);
%
%   Though this looks like 1 complex number, these components of Y should
%   not mix, so the imaginary component will be tacked onto the end of the
%   measurement list to keep them separate.
%
%
%   Example: If data = [1, 10, 100; exp(1), 10*exp(1), 100*exp(1)];
%
%   then LOGMEAN(data) yields [3.6109, 1.3083, -.9943; 3.6109, 1.3083,
%   -.9943].
%
% See Also: LOWPASS, HIGHPASS.
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
dims = size(data_in);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data_in) > 2);
isZ=any(~isreal(data_in(:)));

%% N-D Input.
if NDtf
    data_in = reshape(data_in, [], Nt);
end

%% Perform Logmean.
Phi_0=mean(data_in, 2);
X=bsxfun(@rdivide, data_in, Phi_0);

if ~isZ % All Real
    data_out = -log(X);
else
    Y_Rytov_Re=-log(abs(X));
    Y_Rytov_Im=-angle(X);
    
    data_out=cat(1,Y_Rytov_Re,Y_Rytov_Im);
    dims(1)=2*dims(1);
end


%% Fix any NaNs
data_out(~isfinite(data_out))=0;


%% N-D Output.
if NDtf
    data_out = reshape(data_out,dims);
end



%
