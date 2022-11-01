function [dataR,NRs,beta]=Nuisance_Regress(data,fov,mask,Tregs)
%
% This code calculates regressors as defined by the mean
% time trace witin the intersection of the DOT fov and the given set of
% masks. Generally, these masks delineate the whole brain signal, the
% signals from each hemisphere, or from other regions of interest.
% The mean data time trace within the masks are then simultaneously 
% regressed from the data. If y_{r} is the signal to be
% regressed out and y_{in} is a data time trace, then the output is 
% the least-squares regression: 
%   y_{out} =  y_{in} - y_{r}(<y_{in},y_{r}>/|y_{r}|^2). 
% Additionally, the correlation coefficient is given by:
%   R=(<y_{in},y_{r}>/(|y_{in}|*|y_{r}|)).
% Inputs:   Data-   assumed to be a 4D matrix of data with time as the last
%                   dimension.
%           fov-    a binary mask deliniating the region of data to be
%                   included.
%           mask-   structure containing binary masks of the set of masks
%                   (N=1...Nmask) to be used as the bases of regression.
%           Tregs-  optional input containing time traces to regress. This
%                   must have time as the last dimension, as with the data 
%                   input. The number of elements in the 1st dimension are 
%                   the number of regressors.
% Outputs:  dataR-  the regressed data
%           NRs-    matrix containing nuisance regressors
%           beta-   map of regression coefficient for each regressor with
%                   data

% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht, Zachary E. Markow
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

disp('Calculating regressors')
%% Prepare data
[Nx,Ny,Nz,Nt]=size(data);
data=reshape(data,[],Nt);
if isstruct(mask)
    doSR=1; 
    masks=fieldnames(mask);
    Nmasks=size(masks,1);
    NRs=zeros(Nt,Nmasks);
else
    doSR=0;
    NRs=[];
end

%% Generate Regressors
if doSR                               % Spatial mask-based regressors
for j=1:Nmasks
    foo=mask.(masks{j}).*fov;         % grab mask for regression and limit to fov
    NRs(:,j)=mean(data(foo==1,:),1)'; % Ave time traces within mask
end
end

% Temporal regressors
if exist('Tregs','var'), NRs=cat(1,NRs,Tregs'); end 

%% Pseudoinverse for LS regression
regs_pi=pinv(NRs);
beta=regs_pi*data';

%% Regress and return to full volume
dataR=(data'-NRs*beta)';
dataR=reshape(dataR,Nx,Ny,Nz,Nt);  
beta=reshape(permute(beta,[2,1]),Nx,Ny,Nz,[]);