%% Tikhonov_invert_Amat
% Inverts a sensitivity matrix.
%
%% Description
% |iA = Tikhonov_invert_Amat(A, lambda1, lambda2)| allows the user to
% specify the values of the |lambda1| and |lambda2| parameters in the
% inversion calculation. |lambda2| for spatially-variant regularization,
% is optional.
%
%% See Also
% <smooth_Amat_help.html smooth_Amat> | <reconstruct_img_help.html
% reconstruct_img> | <FindGoodMeas_help.html FindGoodMeas>