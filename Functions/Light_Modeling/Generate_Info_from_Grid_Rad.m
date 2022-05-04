function [info]=Generate_Info_from_Grid_Rad(grid,rad,params,info)
% 
%
% This function generates info structures optodes and pairs from the 
% data in grid and radius files. 
% params can be used to pass in mod type (default is 'CW' but can be the
% modulation frequency if fd) and wavelength(s) of the data in field 'lambda'.
%
% The params input can be used to populate the following information:
%   params.lambda       Wavelengths of the light. Any number of comma-
%                       separated values is allowed. Default: [750,850].
%   params.Mod          Modulation type or frequency. Can be 'CW or 'FD' or
%                       can be the actual modulation frequency (e.g., 0 or
%                       200) in MHz.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
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



%% Parameters and Initialization
if ~exist('info','var'), info=struct;end
if ~exist('params','var'), params=struct;end
if ~isfield(params,'lambda'), params.lambda=[750,850];end
if ~isfield(params,'Mod'), params.Mod='CW';end
Nwl=length(params.lambda);
Nm=size(rad.meas,1);


%% Populate info.optodes structure
if isfield(grid,'name'),info.optodes.CapName=grid.name;end
if isfield(grid,'dpos3'),info.optodes.dpos3=grid.dpos3;end
if isfield(grid,'dpos2')
    info.optodes.dpos2=grid.dpos2;
elseif isfield(grid,'dpos')
    info.optodes.dpos2=grid.dpos;
else
    info.optodes.dpos2=grid.dpos3;
end
if isfield(grid,'spos3'),info.optodes.spos3=grid.spos3;end
if isfield(grid,'spos2')
    info.optodes.spos2=grid.spos2;
elseif isfield(grid,'spos')
    info.optodes.spos2=grid.spos;
else
    info.optodes.spos2=grid.spos3;
end
if ~isfield(info.optodes,'plot3orientation')
    info.optodes.plot3orientation.i='R2L';
    info.optodes.plot3orientation.j='P2A';
    info.optodes.plot3orientation.k='D2V';
end

%% Populate info.pairs table
info.pairs=struct;
info.pairs.Src=repmat(rad.meas(:,1),Nwl,1);
info.pairs.Det=repmat(rad.meas(:,2),Nwl,1);
if isfield(rad,'NN'),info.pairs.NN=repmat(rad.NN,Nwl,1);end
info.pairs.WL = nan(size(info.pairs.Src));
info.pairs.lambda = nan(size(info.pairs.Src));
info.pairs.Mod=repmat(params.Mod,Nm*Nwl,1);
info.pairs.r2d=nan(size(info.pairs.Src));
info.pairs.r3d=nan(size(info.pairs.Src));

for j=1:Nwl
    info.pairs.WL(Nm*(j-1)+1:Nm*j) = ones(Nm,1).*j;
    info.pairs.lambda(Nm*(j-1)+1:Nm*j)=params.lambda(j);
end

for j=1:size(info.pairs.Src,1)
    info.pairs.r2d(j)=norm(info.optodes.spos2(info.pairs.Src(j),:)-...
                info.optodes.dpos2(info.pairs.Det(j),:));
    info.pairs.r3d(j)=norm(info.optodes.spos3(info.pairs.Src(j),:)-...
                info.optodes.dpos3(info.pairs.Det(j),:));
end

%