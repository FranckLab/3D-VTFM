function [s, n, ui] = calculateSurfaceUi(varargin)
%% [s, n, ui] = calculateSurfaceUi(s0, n0, u) is a function that calculates
% surface grid coordinates from the reference surface grid, s0, and
% consecutive finite displacements, u, that were applied by the external
% deformation
%
% INPUTS
% -------------------------------------------------------------------------
%   s0: cell array containing the surface coordinate grid defined by a
%      rectangular meshgrid (can be non-uniform)
%      s0{1}{1} = x-coordinates of nodes
%      s0{1}{2} = y-coordinates of nodes
%      s0{1}{3} = z-coordinates of nodes
%   n0: cell array containing the normal vectors for each point on the 
%      surface
%       n0{1}{1} = x-component of normal vector
%       n0{1}{2} = y-component of normal vector
%       n0{1}{3} = z-component of normal vector
%   u: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix for each
%      time point
%         (components in x,y,z)
%         u{time}{1} = displacement in x-direction
%         u{time}{2} = displacement in y-direction
%         u{time}{3} = displacement in z-direction
%         u{time}{4} = magnitude
% OUTPUTS
% -------------------------------------------------------------------------
%   s: cell array containing the surface coordinate grid defined by a
%      rectangular meshgrid (can be non-uniform) calculated by deforming
%      reference surface, s0, by a known displacement field, u
%      s{time}{1} = x-coordinates of nodes
%      s{time}{2} = y-coordinates of nodes
%      s{time}{3} = z-coordinates of nodes
%   n: cell array containing the normal vectors for each point on the 
%      surface
%       n{time}{1} = x-component of normal vector
%       n{time}{2} = y-component of normal vector
%       n{time}{3} = z-component of normal vector
%   ui: displacement field vector interpolated on the surface grid point with 
%       spacing dm. Format: cell array, each containing a 2D matrix for each
%       time point
%         (components in x,y,z)
%         ui{time}{1} = displacement in x-direction on the surface
%         ui{time}{2} = displacement in y-direction on the surface
%         ui{time}{3} = displacement in z-direction on the surface
%         ui{time}{4} = magnitude
% 
% If used please cite:
% Toyjanova J., Bar-Kochba E., López-Fagundo C., Reichner J., 
% Hoffman-Kim D, Franck, C. (2014) High Resolution, Large Deformation 3D
% Traction Force Microscopy. PLoS ONE 9(4): e90976. 
% doi:10.1371/journal.pone.0090976

[s0,n0,u,dm] = parseInputsUI(varargin{:});

maxTime = length(u);

s = cell(maxTime + 1,1);
ui = cell(maxTime,1);
n = cell(maxTime + 1,1);
s{1} = s0;
n{1} = n0;

for i = 1:maxTime
    
    % interpolate displacements to surface
    ui_ = interpolate2Surface(u(i),dm,s(i),[1 1 1],0);
    
    % check for direction of the normals
    if mean2(n{i}{3}) > 0,
        ui{i} = ui_{1};
    else
        ui{i} = cellfun(@(x) -x, ui_{1},'UniformOutput',0);
        n{i}{3} = -n{i}{3};
    end
    ui{i}{4} = sqrt(ui{i}{1}.^2 + ui{i}{2}.^2 + ui{i}{3}.^2);
    
    % move surface by displacements
    s{i+1} = mapSurface(s{i},ui{i},s0);
    n{i+1} = calculateNormals(s{i+1});
    
end
end

function s = mapSurface(s,ui,s0)

for i = 1:3, s{i} = s{i} - ui{i}; end 
F = TriScatteredInterp(s{1}(:),s{2}(:),s{3}(:));

s{3} = F(s0{1},s0{2});
s{1} = s0{1}; s{2} = s0{2};
s{3} = inpaint_nans(s{3});
end

function [n, tri] = calculateNormals(s)
% tr:  delauney tiranguluzation of the surface. see TriRep()
% s{1}, s{2}: meshgrid that the tractions are defined on
sSize = size(s{1});
dt = delaunay(s{1}(:),s{2}(:));
tri = triangulation(dt,[s{1}(:) s{2}(:) s{3}(:)]);
n_ = vertexNormal(tri);

n = cell(1,3);
for i = 1:3, n{i} = reshape(n_(:,i),sSize); end


end

function varargout = parseInputsUI(varargin)
s0 = varargin{1};
n0 = varargin{2};
u = varargin{3};

if length(s0) == 1, s0 = s0{1}; end
if length(n0) == 1, n0 = n0{1}; end

mSpacing = s0{2}(2) - s0{1}(1);

varargout{      1} = s0;
varargout{end + 1} = n0;
varargout{end + 1} = u;
varargout{end + 1} = mSpacing;

end

function ui = interpolate2Surface(varargin)
% function ui = interpolateDisplacements(u,dm,s,um2vxl,OS)
%               interpolates 3D volumetric displacements on the surface of 
%               the volume. Results are in pixels 

% INPUTS
% -------------------------------------------------------------------------
%   u: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix for each
%      time point (components in x,y,z)
%          u{time}{1} = displacement in x-direction
%          u{time}{2} = displacement in y-direction
%          u{time}{3} = displacement in z-direction
%          u{time}{4} = magnitude
%   dm: the spacing of the meshgrid generated in FIDVC
%   s: cell array containing the surface coordinate grid defined by a
%      rectangular meshgrid (can be non-uniform)
%          s{1}{1} = x-coordinates of nodes
%          s{1}{2} = y-coordinates of nodes
%          s{1}{3} = z-coordinates of nodes
%  um2vxl: micrometer to voxel ratio. Default: [1,1,1]
%  OS: boolean - 0: for reference grid, 1: for deformed grid
%
%
% OUTPUTS
% -------------------------------------------------------------------------
%   ui: displacement field vector interpolated on the surface grid point with 
%       spacing dm. Format: cell array, each containing a 2D matrix for each
%       time point
%         (components in x,y,z)
%         ui{time}{1} = displacement in x-direction on the surface
%         ui{time}{2} = displacement in y-direction on the surface
%         ui{time}{3} = displacement in z-direction on the surface
%         ui{time}{4} = magnitude

[u,m,s,um2vxl,OS] = parseInputs(varargin{:});

method = 'linear';
maxTime = length(u);
ui = cell(maxTime,1);

for i = 1:maxTime
    ui{i} = funInterpolate2Surface(u{i},m,s{i+OS},method,um2vxl);
end


end

function ui = funInterpolate2Surface(u,m,s,method,um2vxl)

if iscell(u), nComponents = length(u); 
else u = {u};
end

nComponents = length(u);
if nComponents == 4, nComponents = 3; end
ui = cell(1,nComponents);
for i = 1:length(ui);
    ui{i} = interp3(m{1},m{2},m{3},u{i},s{1},s{2},s{3},method);
    ui{i} = inpaint_nans(ui{i})*um2vxl(i);
end

switch length(ui)
    case 1, ui = ui{1};
    case 2, ui{3} = sqrt(ui{1}.^2 + ui{2}.^2);
    case 3, ui{4} = sqrt(ui{1}.^2 + ui{2}.^2 + ui{3}.^2);
end

end

function varargout = parseInputs(varargin)
u = varargin{1};
mSpacing = varargin{2};

if iscell(u{1}) == 0
    mSize = size(u{1});
else
    mSize = size(u{1}{1});
end
s = varargin{3};
um2vxl = varargin{4};
OS = varargin{5};


idx = cell(1,3);
for i = 1:3; idx{i} = (1:mSpacing:mSize(i)*mSpacing); end
[m{1}, m{2}, m{3}] = meshgrid(idx{:});

varargout{    1} = u;
varargout{end+1} = m;
varargout{end+1} = s;
varargout{end+1} = um2vxl;
varargout{end+1} = OS;

end
