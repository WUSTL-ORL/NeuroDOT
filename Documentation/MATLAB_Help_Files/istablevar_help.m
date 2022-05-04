%% istablevar
% Checks if a variable or list of variables exist in a table.
% 
%% Description
% |istablevar(table, varname)| performs the same logical check on tables as
% <matlab:doc('isfield') |isfield|> does for structure arrays, by checking
% if the variable(s) |varname| exist within the table |table|. |varname|
% must be either a string or a 1D cell array of strings. A horizontal cell
% array will give a logical array for each cell array element, while a
% vertical array will be equivalent to |all(istablevar)|.
% 
% Example: If |info.pairs| is a table with variable names |Src|, |Det|,
% |WL|, |NN|, and |r|, then,
istablevar(info.pairs, {'WL'})
%%
% gives
% 
%   ans = 
%           1
% 
% and with a horizontal cell array,
istablevar(info.pairs, {'WL', 'NN'})
%%
% gives
% 
%   ans = 
%           1   1
% 
% However, note that the component function |strcmp| in |istablevar| is case-sensitive, so,
% 
istablevar(info.pairs, {'wL', 'Nn'})
%% 
% will give 
%
%   ans =
%           0   0
% 
% and with a vertical input array will result in a single output, only returning |true| if both strings are matches.
% 
istablevar(info.pairs, {'WL'; 'NN'})
%% 
%   ans =
%           1
% 
%% See Also
% <matlab:doc('is*') is*> | <matlab:doc('isfield') isfield> |
% <matlab:doc('all') all>