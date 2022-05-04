function stats=nmi(A,B)
%
% This function calculates various information-theoretical metrics for the
% passed variables including: H-entropy, Mab-mutual information,
% NMI-normalized mutual information.
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

%% 
N=numel(A);
A=A(:);
B=B(:);

l=min(min(A),min(B));
A=A-l+1;
B=B-l+1;
u=max(max(A),max(B));

idx=[1:N]';
Ma=sparse(idx,A,1,N,u,N);
Mb=sparse(idx,B,1,N,u,N);


%% Prob Distributions
Pa=mean(Ma,1);          % Individual distributions
Pb=mean(Mb,1);
Pab=nonzeros(Ma'*Mb/N); % Joint distribution


%% Entropies
Ha=-dot(Pa,log2(Pa+eps));
Hb=-dot(Pb,log2(Pb+eps));
Hab=-dot(Pab,log2(Pab+eps));


%% Mutual Information
Mab=Ha+Hb-Hab;
NMI=sqrt((Mab/Ha)*(Mab/Hb));

%% outputs
stats.NMI=NMI;
stats.Mab=Mab;
stats.Ha=Ha;
stats.Hb=Hb;
stats.Pa=Pa;
stats.Pb=Pb;