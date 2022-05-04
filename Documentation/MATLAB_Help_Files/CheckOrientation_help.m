%% CheckOrientation
% Interprets 3D orientation information.
%
%% Description
% |[xdir, ydir, zdir] = CheckOrientation(info)| checks the directions of the
% data encoding in |info.optodes.plot3orientation| to make it compatible
% with 3D cap visualizations. The default reverses the x-axis and leaves
% the other two alone. The directions are returned in |xdir|, |ydir|, and
% |zdir|.
% 
%% See Also
% <PlotCapData_help.html PlotCapData>