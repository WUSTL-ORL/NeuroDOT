function [pFDR,Htest]=FDR(Pin,alpha,mask)
%
% This function runs a False Discovery Rate (FDR) test on a set of N
% p-values. Per Benjamini & Hochberg 1995, the FDR procedure is:
%   (0) Obtain a set of p-values for a test
%   (1) Sort the p-values from lowest-to-highest
%   (2) Calculate pFDR(i,alpha) = i*alpha/N; where is is the sorted index
%   (3) Test that Pin<=pFDR. If so, Htest=1.
% Pin is taken only in arrary elements where mask==1. Pin is assumed to be
% in the range [0,1];


%% Paramters
if ~exist('alpha','var'), alpha=0.05;end
if ~exist('mask','var'), mask=ones(size(Pin),'single');end
% NZ=intersect(find(Pin<1),find(Pin>0));
NZ=find(mask==1);
Nps=size(NZ,1);
Pvals=Pin(NZ);
pFDR=zeros(size(Pin),'single');
Htest=zeros(size(Pin),'single');


%% Sort p-values
[~,idx]=sort(Pvals);

%% Calc pFDR and test
pFDR(NZ(idx))=([1:Nps]).*(alpha./Nps);
Htest(NZ(idx))=pFDR(NZ(idx))>=Pin(NZ(idx));
