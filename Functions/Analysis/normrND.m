function data=normrND(data)
%
% This function returns a row-normed matrix. It is assumed that the matrix
% is 2D. Updated for broader compatability.

%%
% data=bsxfun(@rdivide,data,vecnorm(data,2,2));
dataNorm=sqrt(sum(data.^2,2));
data=bsxfun(@rdivide,data,dataNorm);
data(~isfinite(data))=0;