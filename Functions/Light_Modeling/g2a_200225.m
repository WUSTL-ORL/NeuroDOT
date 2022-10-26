function [A,Gsd]=g2a_200225(Gs,Gd,dc,dim,flags)

% g2a(): Make A-Matrix with Adjoint Method and Rytov Approximation
% Gs    Green's functions for sources
% Gd    Green's functions for detectors
% dc    Diffusion coefficient
% dim   space meta data structure
% flags structure containing parameters, including Hz for mod freq of
%           source light.
% 
% %%200225 Added functionality to allow A for a subset of measurements as
%           listed in measurement list in flags.infoA.pairs.
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


%% Parameters and Initialization
if ~isfield(flags,'Hz'),flags.Hz=0;end
if ~isfield(flags,'compute_mua'), flags.compute_mua=1;end
if ~isfield(flags,'compute_mus'), flags.compute_mus=0;end
if ~isfield(flags,'infoA'),disp(['Must include measurement list']);return; end
    
if flags.Hz~=0 && flags.compute_mua && flags.compute_mus
    % If both, make A as a structure
    flagsD=flags;
    flagsD.compute_mua=0;
    [A.mus.A,A.mus.Gsd]=g2a_200225(Gs,Gd,dc,dim,flagsD); % Scattering
    
    flagsA=flags;
    flagsA.compute_mus=0;
    [A.mua.A,A.mua.Gsd]=g2a_200225(Gs,Gd,dc,dim,flagsA); % Absorption
    Gsd=[];
    return
end
info=flags.infoA; % Herein is the measurement list
numpt=size(Gd,3); % Gx are lambda-SD-vox
measnum=size(info.pairs,1);
disp(['Initializing A: ',num2str(measnum),' Measurements by ',...
    num2str(numpt),' voxels'])
A=zeros(measnum,numpt,'single');
Gsd=zeros(measnum,1);


%% Calculate Gsd
% Acw_denominator = Gs(d)
disp('Calculating Gsd')
tic
[~,Gdidx]=max(abs(Gd),[],3); %Gs(d)
for m=1:measnum
    Gsd(m)=Gs(info.pairs.WL(m),info.pairs.Src(m),Gdidx(info.pairs.WL(m),...
        info.pairs.Det(m)));
end
toc


%% Adjoint Formulation and Normalization: 
if flags.Hz==0 || flags.compute_mua == 1    
    disp(['Creating Adjoint formulation for mua'])
    tic
    for m=1:measnum
        A(m,:)=squeeze(Gs(info.pairs.WL(m),info.pairs.Src(m),:)).*...
            squeeze(Gd(info.pairs.WL(m),info.pairs.Det(m),:)).*...
            (dim.sV^3./dc(info.pairs.WL(m),:))./Gsd(m);
    end
    toc

elseif flags.compute_mus              
        disp(['Creating Adjoint formulation for musp'])
        tic
    for m=1:measnum
        foos = zeros(dim.nVx,dim.nVy,dim.nVz); % careful: may be large.
        food = foos;
        foos(dim.Good_Vox) = squeeze(Gs(info.pairs.WL(m),info.pairs.Src(m),:));
        [gGsx,gGsy,gGsz] = gradient(foos,dim.sV);
        
        food(dim.Good_Vox) = squeeze(Gd(info.pairs.WL(m),info.pairs.Det(m),:));
        [gGdx,gGdy,gGdz] = gradient(food,dim.sV);
        
        A(m,:) =-(gGsx(dim.Good_Vox).*gGdx(dim.Good_Vox)+...
                    gGsy(dim.Good_Vox).*gGdy(dim.Good_Vox)+...
                    gGsz(dim.Good_Vox).*gGdz(dim.Good_Vox)).*...
                    (dim.sV^3./dc(info.pairs.WL(m),:))./Gsd(m);
    end   
    toc
else A=[];Gsd=[];return
end


%% Normalize Rytov with Gsd and interp with vol and diff coef
% disp(['Normalizing Rytov'])
% tic
% A=bsxfun(@rdivide,A,Gsd);
% toc
%     
% disp(['>Normalizing for Discretized Space'])
% tic
% for m=1:measnum
%     f=dim.sV^3./dc(info.pairs.WL(m),:);
%     A(m,:)=f.*squeeze(A(m,:));
% end
% toc

%% Remove NaN's
A(isnan(A))=0;