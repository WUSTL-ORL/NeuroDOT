function Par_tt = Resample_Vol_to_Parcel(volume, Parcels, meshL, meshR, dim)

% Resample_Vol_to_Parcel Resamples data to a single time trace for each parcel in the
% field of view of a specified volume.
% 
%   Resample_Vol_to_Parcel(volume, Parcels, meshL, meshR, dim)
%
%   Inputs:
%
%       volume: Volumetric data 
%           - Must be 4-D, containing time information
%
%       Parcels: The structure containing the parcel data
%           - Must contain "CtxL" and "CtxR" fields
%
%       meshL: Left Hemisphere cortical mesh
%       
%       meshR: Right Hemisphere cortical mesh
%       
%       dim: Info associated with "volume" 
%
%   Outputs:
%
%       Par_tt: Array of # parcels x # of time points. One time trace per
%       parcel.
%
% Dependencies: VOL2SURF_MESH
% See Also: SCRIPT_FOR_VIEWING_PARCELS
%
% Copyright (c) 2024 Washington University 
% Created  by: Emma Speh and Dalin Yang
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
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


%% Interpolate volumetric data onto the surface mesh 
mapL = vol2surf_mesh(meshL, volume, dim); 
mapR = vol2surf_mesh(meshR, volume, dim); 
ParNum=max(cat(1,Parcels.CtxL,Parcels.CtxR)); 
% Initialize time traces per parcel
Par_tt = zeros(length(unique(cat(1,Parcels.CtxL,Parcels.CtxR))), size(volume,4));   


%% Average volumetric data within each parcel
indexpar = 0;
for Par=1:ParNum
    idxP = find(Parcels.CtxL==Par);       % find the surface nodes within each parcel
    idx1 = any(mapL.data,2);              % find the surface nodes within the field of view
    idxD = find(idx1);
    idx = intersect(idxP,idxD);           % intersect parcel data and parcel field of view
    if (isempty(idx)==0)
        indexpar = indexpar+1;
        mu=mean(mapL.data(idx,:));        % average the volumetric data within each parcel    
        Par_tt(indexpar,:) = mu;
    end
    idxP = find(Parcels.CtxR==Par);
    idx1 = any(mapR.data,2);              % find the surface nodes within the field of view
    idxD = find(idx1);
    idx = intersect(idxP,idxD); 
    if (isempty(idx)==0)
        indexpar = indexpar+1;
        mu=mean(mapR.data(idx,:));        % average the volumetric data within each parcel  
        Par_tt(indexpar,:) = mu;
    end
end
deleteRows = all(Par_tt==0,2);
Par_tt(deleteRows,:)=[];

end


%%

