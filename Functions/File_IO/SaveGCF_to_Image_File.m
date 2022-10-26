function SaveGCF_to_Image_File(fn, kill, fig_handle, format_type, res_dpi, use_CMYK)
%
%% Program Description
%
% SaveGCF_to_Image_File(fn, kill, fig_handle, format_type, res_dpi, use_CMYK)
%
% This function outputs a figure to the filename fn in the specified
% format_type (default is 'png'). See below for a list of common acceptable
% options for format_type.  For a more complete list of these options, see
% the documentation on MATLAB's 'print' function.  After saving, the figure
% is closed if kill is unspecified or set to 1 (default); set kill to 0 to
% avoid closing the figure.  res_dpi is the image resolution in dots per
% inch; 600 dpi is used by default, and 0 dpi can be specified to use the
% same resolution as the screen resolution.  use_CMYK is true by default,
% but according to MATLAB's 'print' function documentation, the
% corresponding '-cmyk' flag only changes the color scheme CMYK for EPS and
% PS files. User can set use_CMYK to false if RGB is desired instead of
% CMYK in an EPS or PS file.  Set any input variable except fn to the empty
% array [] to use the default setting for that input variable.
%
% Based on SaveGCF_to_PNG created by Adam T. Eggebrecht, extended into this
% function by Zachary E. Markow to provide additional formatting options.
%
% List of common options for format_type recognized by this function:
%
% - 'png', 'PNG', or '-dpng' (PNG file)
% - 'eps', 'EPS', or '-depsc' (EPS file with Level 3 color)
% - 'jpg', 'JPG', 'jpeg', 'JPEG', or '-djpeg' (JPEG file)
% - 'tiff', 'TIFF', or '-dtiff' (compressed TIFF file)
% - '-dtiffn' (uncompressed TIFF file)
% - See MATLAB's 'print' function documentation for more options.
%


%% Set defaults.

if (~exist('kill','var')) || isempty(kill)
    kill=1;
end
if (~exist('fig_handle','var')) || isempty(fig_handle)
    fig_handle=gcf;
end
if (~exist('format_type','var')) || isempty(format_type)
    format_type = 'png';
end
if (~exist('res_dpi','var')) || isempty(res_dpi)
    res_dpi = 600;
end
if (~exist('use_CMYK','var')) || isempty(use_CMYK)
    use_CMYK = true;
end


%% Convert shortcut versions of format_type if necessary.

format_type_lower = lower(format_type);

switch format_type_lower
    case 'png'
        format_type_actual = '-dpng';
    case 'eps'
        format_type_actual = '-depsc';
    case {'jpg','jpeg'}
        format_type_actual = '-djpeg';
    case 'tiff'
        format_type_actual = '-dtiff';
    otherwise
        format_type_actual = format_type;
end


%% Save file.

res_str = ['-r' num2str(round(res_dpi))];

pause(1)
set(fig_handle, 'InvertHardCopy', 'off');
set(fig_handle,'PaperPositionMode','auto')
if use_CMYK
    print(fn, format_type_actual, res_str, '-cmyk')
else
    print(fn, format_type_actual, res_str)
end


%% Close figure if user desires.

if kill, close(fig_handle); end
pause(1)


end  % End of function.