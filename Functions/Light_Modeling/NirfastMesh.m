function nirfast_mesh=NirfastMesh(mask,meshname,param,Mode)
%
% Mask is the segmented image.
% Mode allows for faster mesh generation for models that will not be used
% for light modeling.
% Key params from cGal help for surface facets and mesh cells:
% from http://doc.cgal.org/latest/Mesh_3/index.html#Chapter_3D_Mesh_Generation
%   facet_angle: This parameter controls the shape of surface facets. 
%                   Actually, it is a lower bound for the angle 
%                   (in degree) of surface facets. When boundary 
%                   surfaces are smooth, the termination of the 
%                   meshing process is guaranteed if the angular 
%                   bound is at most 30 degrees.
%                   Recommend: 25 to 30, not more.
%   facet_size:  This parameter controls the size of surface facets. 
%                   Actually, each surface facet has a surface Delaunay 
%                   ball which is a ball circumscribing the surface 
%                   facet and centered on the surface patch. The 
%                   parameter facet_size is either a constant or a 
%                   spatially variable scalar field, providing an 
%                   upper bound for the radii of surface Delaunay balls.
%                   Recommended ~ 0.15. SMALLER = MORE DENSE
%   facet_distance:    This parameter controls the approximation error of 
%                   boundary and subdivision surfaces. Actually, it 
%                   is either a constant or a spatially variable scalar 
%                   field. It provides an upper bound for the distance 
%                   between the circumcenter of a surface facet and the 
%                   center of a surface Delaunay ball of this facet.
%                   Recommended ~ 0.05. Smaller values will more precisely
%                   match the input surfaces.
%   facet_topology:    This parameters controls the set of topological 
%                   constraints which have to be verified by each 
%                   surface facet. By default, each vertex of a surface 
%                   facet has to be located on a surface patch, on a 
%                   curve segment, or on a corner. It can also be set 
%                   to check whether the three vertices of a surface 
%                   facet belongs to the same surface patch. This has 
%                   to be done cautiously, as such a criteria needs 
%                   that each surface patches intersection is an input 
%                   1-dimensional feature.
%   cell_radius_edge_ratio:  This parameter controls the shape of mesh 
%                   cells (but can't filter slivers, as we discussed 
%                   earlier). Actually, it is an upper bound for the 
%                   ratio between the circumradius of a mesh 
%                   tetrahedron and its shortest edge. There is a 
%                   theoretical bound for this parameter: the Delaunay 
%                   refinement process is guaranteed to terminate for 
%                   values of cell_radius_edge_ratio bigger than 2.
%                   The mesh quality referred to here is measured by 
%                   the radius edge ratio of surface facets end mesh 
%                   cells, where the radius edge ratio of a simplex 
%                   (triangle or tetrahedron) is the the ratio between 
%                   its circumradius and its shortest edge.
%                   Recommended >2, and at or under 4
%   cell_size:   This parameter controls the size of mesh tetrahedra. 
%                   It is either a scalar or a spatially variable scalar 
%                   field. It provides an upper bound on the circumradii 
%                   of the mesh tetrahedra.
%                   Recomended ~ 0.2. SMALLER = MORE DENSE
%
%
%
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
if ~exist('param','var'),param='hd';end
if ~isstruct(param),type=param;param=struct;end
if ~exist('type','var'),type='hd';end
if ~ischar(type), value=type;type='custom';end
if ~exist('Mode','var'),Mode=0;end
if ~isfield(param,'CheckMeshQuality'),param.CheckMeshQuality=0;end

stack = uint8(mask);

% Pixel size and slice spacing in mm 
if isfield(param,'info')
    sx=param.info.mmx;
    sy=param.info.mmy;
    sz=param.info.mmz;
    if ~isfield(param,'Offset')
    param.Offset = param.info.center.*((-1).^((param.info.mmppix>0)));
    end
else
    sx = 1; sy = 1; sz = 1;
    param.Offset = [0 0 0];
end
param.medfilter = 0; 
param.PixelSpacing(1) = sx;
param.PixelSpacing(2) = sy;
param.SliceThickness  = sz;

% For UHD (assume input is [1,1,1]), try:
% param.cell_size=0.5;
% param.facet_distance=0.5;

% Specify global tetrahedron size
if ~isfield(param,'cell_size')
    switch type
        case 'hd'
            param.cell_size = 2.0 * min(sx, min(sy,sz));
        case 'ld'
            param.cell_size = 5.5 * min(sx, min(sy,sz));
        case 'custom'
            param.cell_size = value * min(sx, min(sy,sz));
    end
end

if ~isfield(param,'facet_size')
param.facet_size = param.cell_size; % for surf and vol to match
end
    
if ~isfield(param,'cell_radius_edge')
param.cell_radius_edge = 3.0; % keep between [2,4]
end

if ~isfield(param,'facet_distance')
param.facet_distance = 3.0; 
end

if ~isfield(param,'facet_angle')
param.facet_angle = 25.0;  
end

param.special_subdomain_label = 0;
param.special_subdomain_size  = 0;
param.tmppath = cd;

%% Run mesh generator
disp('<<<Generating initial mesh:')
disp([param])

t=tic;
[e,p] = RunCGALMeshGenerator(stack,param);
disp(['>>Mesh Generation took ',num2str(toc(t)),' seconds'])


%% Fix Offset if non-zero and x-dim is flipped
if any(param.Offset)
    if length(unique(sign(param.Offset)))>1
        if param.Offset(1)>0
            p(:,1)=-p(:,1)+2*param.Offset(1);
        end
    end
end


%% Call conversion to nirfast mesh
if Mode
    nirfast_mesh.elements=double(e(:,1:4));
    nirfast_mesh.nodes=p;
else
    mesh.ele = e;
    mesh.node = p;
    mesh.nnpe = 4;
    mesh.dim = 3;
    
    % Change mesh type based on type of analysis being performed. Accepted
    % values are: 'stnd', 'spec', 'stnd_spn', 'fluor'
    disp('Generating Nirfast-compliant mesh')
    if ~isfield(param,'meshtype')
        meshtype = 'stnd';
    else
        meshtype=param.meshtype;
    end
    solidmesh2nirfast(mesh,[meshname],meshtype);
    
    %% Load nirfast mesh
    nirfast_mesh = load_mesh([meshname]);
end

%% Check Mesh Quality and visualize
if param.CheckMeshQuality
disp(['< Assessing mesh quality']);
[vol,q,q_area,status]=CheckMesh3D(nirfast_mesh.elements,nirfast_mesh.nodes);

vol=abs(vol);
vm=min(vol(:));
vM=max(vol(:));
dv=(vM-vm)/100;
qm=min(q(:));
qM=max(q(:));
dq=(qM-qm)/100;
figure('Color','w');
subplot(1,2,1);histogram(vol(:),[(vm-dv):dv:(vM+dv)]);
xlabel('Volume');ylabel('N elements');
subplot(1,2,2);histogram(q(:),[(qm-dq):dq:(qM+dq)]);
xlabel('Simplex Quality Ratio');ylabel('N elements');
switch status.solid
    case 0
        disp(['< Mesh quality: acceptable']);
    case 2
        disp(['< Mesh quality: some elements have close to zero volume']);
    case 4
        disp(['< Mesh quality: some elements have very low quality']);
    case 8
        disp(['< Mesh quality: serious tetrahedron connectivity issue,',...
            'some of faces are shared by more than two tetrahedrons']);
end
end