function [fwhm,stats]=Find_FWHM(Vin,prior,th,type)
%
% This function calculates the full width half maximum of an activation
% given a prior. If desired, a threshold can be defined (default = 0.5).
% Inputs:
%   Vin     Volumetric data
%   prior   Matlab indices of spatial prior to initialize FWHM search
%   th      Relative threshold multiplier (0.5 = 50%)
%   type    keep island 0=with prior, 1=largest, 2=strongest continuous region
%               If 1/2, prior may not lie within returned island.
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

if ~exist('th','var'), th=0.5;end

if exist('prior','var') && ~isempty(prior)
    pVal=Vin(prior(1),prior(2),prior(3));
    pIdx= sub2ind(size(Vin),prior(1),prior(2),prior(3));
else
    [pVal,pIdx]=max(Vin(:));
end
if ~exist('type','var'), type=0;end % 0: with prior, 1: biggest; 2: strongest


%% Find islands
test0=+(Vin>(th*pVal)); % Threshold and Binarize
test1=bwconncomp(test0); % find islands
if test1.NumObjects>1
    keep=zeros(test1.NumObjects,3);
    for j=1:test1.NumObjects
        keep(j,1)=sum(ismember(test1.PixelIdxList{j},pIdx));    % w prior
        keep(j,2)=length(test1.PixelIdxList{j});                % size
        keep(j,3)=sum(Vin(test1.PixelIdxList{j}));              % strength
    end
    test2=zeros(size(test0));
    switch type
        case 0
            test2(test1.PixelIdxList{keep(:,1)==1})=1; % Island w prior
            fwhm=test2.*Vin;
        case 1
            test2(test1.PixelIdxList{keep(:,2)==max(keep(:,2))})=1; % Largest
            test2=test2.*Vin; % Set new prior based on island max
            fwhm=Find_FWHM(test2,[],th,type);
        case 2
            test2(test1.PixelIdxList{keep(:,3)==max(keep(:,3))})=1; % Strongest
            test2=test2.*Vin; % Set new prior
            fwhm=Find_FWHM(test2,[],th,type);
    end
else
    fwhm=zeros(size(test0));
    fwhm(test1.PixelIdxList{1})=Vin(test1.PixelIdxList{1});
end


%% Centroid and peak locations
keep=find(fwhm==max(fwhm(:)));
[x,y,z]=ind2sub(size(Vin),keep);
stats.peak=[x,y,z]; % peak

keep=find(fwhm~=0);
fu=Vin(keep);
[x,y,z]=ind2sub(size(Vin),keep);
r=[x,y,z];
stats.fwhm=max(max(pdist2(r,r)));              % fwhm
stats.Cent=[mean((x-mean(x)).*fu)+mean(x),... % Centroid
            mean((y-mean(y)).*fu)+mean(y),...
            mean((z-mean(z)).*fu)+mean(z)];

stats.fvhm=length(keep); % fvhm