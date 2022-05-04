function data=normcND(data)
%
% This function returns a column-normed matrix. It is assumed that the 
% matrix is 2D.

%%
data=bsxfun(@rdivide,data,vecnorm(data,2,1));