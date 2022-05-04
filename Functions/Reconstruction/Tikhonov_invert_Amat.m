function iA = Tikhonov_invert_Amat(A, lambda1, lambda2)

% TIKHONOV_INVERT_AMAT Inverts a sensitivity matrix.
%
%   iA = TIKHONOV_INVERT_AMAT(A, lambda1, lambda2) allows the user to
%   specify the values of the "lambda1" and "lambda2" parameters in the
%   inversion calculation. lambda2 for spatially-variant regularization,
%   is optional.
%
% See Also: SMOOTH_AMAT, RECONSTRUCT_IMG, FINDGOODMEAS.
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
if ~isreal(A)                   % If complex A, sep into [Re;Im] first
    A=cat(1,real(A),imag(A));
end
[Nm, Nvox] = size(A);
if exist('lambda2','var')
    svr=1;
else
    svr=0;
end


%% Spatially variant regularization
if svr
    ll_0=sum(A.^2,1); % diag(L) = diag(A'A) (shortcut)
    ll=sqrt(ll_0+lambda2*max(ll_0)); % Adjust with Lambda2 cut-off value
    A=bsxfun(@rdivide,A,ll);
end

%% Take the pseudo-inverse.
if Nvox < Nm
    Att = zeros(Nvox, 'single'); % Preallocate
    Att = single(A' * A);
    ss = normest(Att); % normest used because Nvox x Nvox "Att" array is very large.
    penalty = sqrt(ss) .* lambda1;
    iA = (Att + penalty .^ 2 .* eye(Nvox, 'single')) \ A';
else
    Att = zeros(Nm, 'single'); % Preallocate
    Att = single(A * A');
    ss = norm(Att);
    penalty = sqrt(ss) .* lambda1;
    iA = A' / (Att + penalty .^ 2 .* eye(Nm, 'single'));
end


clear Att A ss

%% Undo spatially variant regularization 
if svr
    iA=bsxfun(@rdivide,iA,ll');
end


%
