function SaveGCF_to_PNG(fn,kill,fig_handle)
%
% This function outputs a figure to the filename fn as a png, and 
% closes it

if ~exist('kill','var'),kill=1;end
if ~exist('fig_handle','var'),fig_handle=gcf;end


pause(1)
set(fig_handle, 'InvertHardCopy', 'off');
set(fig_handle,'PaperPositionMode','auto')
print(fn,'-dpng','-r600','-cmyk')

if kill, close; end
pause(1)