function [t,p,dof,cohensD]=TtestWelch(X1,X2)
%
% This function calculates the t-test between 2 sets of data using the
% Welch method that does not assume equal mean or variance or samples. The
% function returns the t-statistic, p-value, and degrees of freedom.
% The input data, X1 and X2, are assumed to be multi-dimensional arrays of
% equal size except for the last dimension that represents
% subjects/instances/etc.
% Cohen's d is also calculated and returned in the 4th output.

%% Prepare data
Xsize1=size(X1);if Xsize1(end)==1, Xsize1(end)=[];end
Xsize2=size(X2);if Xsize2(end)==1, Xsize2(end)=[];end
Ndim1=length(Xsize1);
Ndim2=length(Xsize2);
if Ndim1~=Ndim2
    disp('Data arrays must have same number of dimensions');
    t=zeros(size(min([Ndim1,Ndim2])),'single');
    p=ones(size(min([Ndim1,Ndim2])),'single');
    return
else Ndim=Ndim1;
end

%% Means
mu1=mean(X1,Ndim);
mu2=mean(X2,Ndim);

%% Ns
n1=Xsize1(Ndim);
n2=Xsize2(Ndim);

%% SE2
se1=var(X1,0,Ndim)./n1;
se2=var(X2,0,Ndim)./n2;

%% Unbiased estimator of variance; Satterthwaite
d=sqrt(se1+se2);

%% T-statistic
t=(mu1-mu2)./d;
t(~isfinite(t))=0;

%% Degrees-of-freedom; Welch-Satterthwaite
dof=((se1+se2).^2)./(((se1.^2)./(n1-1))+((se2.^2)./(n2-1)));
dof(~isfinite(dof))=0;

%% P-values: 2-tailed p-value. to 1st approx, halve if want 1-tailed.
p=1-tcdf(abs(t),dof); 
p(~isfinite(p))=0;

%% Cohen's d
cohensD=(mu1-mu2)./sqrt((((n1-1)*n1*se1)+((n2-1)*n2*se2))/(n1+n2-2));
cohensD(~isfinite(cohensD))=0;