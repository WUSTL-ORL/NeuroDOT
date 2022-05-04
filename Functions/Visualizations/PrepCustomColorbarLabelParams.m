function params = PrepCustomColorbarLabelParams(params)
%
% PrepColorbarLabelParams: This function prepares a params structure for
% custom colorbar labeling.  Default values are added for any related
% settings that are missing from params.  However, the inputted params data
% structure must have its Scale field defined (params.Scale must exist).
% If params.Scale does not exist, it will be set to 1.
%
% This function is useful for setting params.cbticks to its default or
% according to custom tick values specified in params.cblabels.  This is
% because for many NeuroDOT plots, the colors of plotted data items (e.g.,
% dots, lines) are determined by applycmap(...) and are not automatically
% tied to a colorbar by MATLAB.  Therefore when "colorbar()" is run to
% display a colorbar, MATLAB's default colorbar is generated, which behaves
% as if displaying data between 0 and 1.  Therefore, to place a tick mark
% halfway up the colorbar (regardless of what data value actually
% corresponds to this halfway mark), we have to tell MATLAB to place a tick
% mark at 0.5, so 0.5 would need to be an entry in the vector
% params.cbticks.
%
% See Also: PLOTSLICES, PLOTINTERPSURFMESH, APPLYCMAP.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht, Zachary E. Markow
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


%% Define some helpful intermediate variables based on inputted params.
if ~isfield(params, 'Scale'), params.Scale=1;end
c_max = params.Scale;
if isfield(params, 'PD')  &&  ~isempty(params.PD)  &&  params.PD
    c_mid = c_max / 2;
    c_min = 0;
else
    c_mid = 0;
    c_min = -c_max;
end


%% If cbticks_in_data_units exists, convert cbticks_in_data_units to cbticks.  Also prepare cblabels accordingly if not already done.

if isfield(params,'cbticks_in_data_units') && (~isempty(params.cbticks_in_data_units))
    params.cbticks = (params.cbticks_in_data_units - c_min) ./ (c_max - c_min);
    if (~isfield(params,'cblabels')) || isempty(params.cblabels)
        params.cblabels = sprintf('%.3g ', params.cbticks_in_data_units);
        params.cblabels = strsplit(params.cblabels);
        params.cblabels = params.cblabels(1:(end-1));  % Remove last '' array element from extra space at end of sprintf result.
    end
end


%% Interpret and (if necessary) create or modify params.cbticks and params.cblabels.

if (~isfield(params, 'cbticks')  ||  isempty(params.cbticks))  &&...
        (~isfield(params, 'cblabels')  ||  isempty(params.cblabels))
    % If both are empty/not here, fill in defaults.
    params.cbticks = [0, 0.5, 1];
    params.cblabels = strtrim(cellstr(num2str([c_min, c_mid, c_max]', 3)));
elseif ~isfield(params, 'cbticks')  ||  isempty(params.cbticks)
    % If only ticks missing...
    if numel(params.cblabels) > 2
        if isnumeric(params.cblabels)
            % If we have numbers, we can sort then scale them.
            params.cblabels = sort(params.cblabels);
            params.cbticks = (params.cblabels - params.cblabels(1))...
                / (params.cblabels(end) - params.cblabels(1));
        else
            params.cbticks = 0:1/(numel(params.cblabels) - 1):1;
        end
    elseif numel(params.cblabels) == 2
        params.cbticks = [0, 1];
    else
        error('*** Need 2 or more colorbar ticks. ***')
    end
elseif ~isfield(params, 'cblabels')  ||  isempty(params.cblabels)
    if numel(params.cbticks) > 2
        % If only labels missing, scale labels to tick spacing.
        scaled_ticks = (params.cbticks - params.cbticks(1)) /...
            (params.cbticks(end) - params.cbticks(1));
        params.cblabels = strtrim(cellstr(num2str([scaled_ticks * (c_max - c_min) + c_min]', 3)));
    elseif numel(params.cbticks) == 2
        params.cblabels = strtrim(cellstr(num2str([c_min, c_max]', 3)));
    else
        error('*** Need 2 or more colorbar labels. ***')
    end
elseif numel(params.cbticks) == numel(params.cblabels)
    % As long as they match in size, continue on.
else
    error('*** params.cbticks and params.cblabels do not match. ***')
end


end  % End of function.