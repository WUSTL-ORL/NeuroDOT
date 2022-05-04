%% Cap_Fitter
% An interactive tool for fitting optode grids onto the skull.
%
%% Description
% |Cap_Fitter| is an interactive GUI that allows the user to load skull
% mesh and cap grid files, re-orient them with an affine transformation,
% and then run a spring-fit algorithm on them. The final result will then
% be saved in a file containing the original "info" structure, the
% resulting one, the mesh they were fit onto, and an "ops" structure
% recording the operations performed.
%
%% Dependencies
% <PlotMeshSurface_help.html PlotMeshSurface> | <PlotCap_help.html PlotCap>
% | <springfit_optodes_help.html springfit_optodes> |
% <operate_optodes_help.html operate_optodes> | <scale_cap_help.html
% scale_cap> | <rotate_cap_help.html rotate_cap> |
% <rotation_matrix_help.html rotation_matrix>