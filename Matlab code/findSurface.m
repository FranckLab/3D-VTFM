function [s,n,angle001] = findSurface(varargin)
% function [s,n,angle001] = calculateSurface(I0 or fileName, dataChannel, dm,
% smoothness, blockSize,plot01) is used find and generate the surface grid
% and normals based on the volumetric image, I0
%
%
% INPUTS
% -------------------------------------------------------------------------
%    I0 or fileName: volumetric image of the first time stack
%    dataChannel: If the volume is a cell array - pick the value of
%                 the volume of beads
%    dm: spacing of the meshgrid
%    smoothness: scalar or vector of length 2 - determines the
%                eventual smoothness of the estimated surface. A larger
%                value here means the surface will be smoother. Smoothness
%                must be a non-negative real number.
%    blockSize: size of the block for averaging the bead intensity (depends
%               on the bead density in the volume)
%    plot01: boolean - choice 0 or 1 for plotting the surface with 
%            surface normals
%
% OUTPUTS
% -------------------------------------------------------------------------
%   s: cell array containing the surface coordinate grid defined by a
%      rectangular meshgrid (can be non-uniform)
%      s{time}{1} = x-coordinates of nodes
%      s{time}{2} = y-coordinates of nodes
%      s{time}{3} = z-coordinates of nodes
%   n: cell array containing the normal vectors for each point on the 
%      surface
%       n{time}{1} = x-component of normal vector
%       n{time}{2} = y-component of normal vector
%       n{time}{3} = z-component of normal vector
% angle001: angle deviation from [0,0,1] surface normal
% 
% If used please cite:
% Toyjanova J., Bar-Kochba E., López-Fagundo C., Reichner J., 
% Hoffman-Kim D, Franck, C. (2014) High Resolution, Large Deformation 3D
% Traction Force Microscopy. PLoS ONE 9(4): e90976. 
% doi:10.1371/journal.pone.0090976

[I, dm, weights, smoothness,blockSize, plot01] = parseInputs(varargin{:});

maxTime = length(I);
s = cell(maxTime,1);
n = cell(maxTime,1);
angle001 = cell(maxTime,1);

for i = 1:maxTime
    [s{i}, centroid] = funCalculateSurface(I{i}, weights, dm, smoothness,blockSize,flag01);
    [n{i}, tri] = funTriInterpNormals(s{i},centroid);
    angle001{i} = calculateNormalAngle(n{i},[0 0 1]);
    if plot01, plotSurface(tri,s{i},n{i},angle001{i}); end
end

end

function [s, centroid] = funCalculateSurface(I0, weights, dm, smoothness, blockSize)
sizeI = size(I0);

% Interpolate surface onto rectangular meshgrid
s1Idx = cell(1,3); s4Idx = cell(1,3);
for i = 1:3,
    s4Idx{i} = (1:dm:(sizeI(i))+1);
    s1Idx{i} = (1:1:(sizeI(i)));
end

[s4{1}, s4{2}] = meshgrid(s4Idx{1},s4Idx{2});

% Threshold image
I = double(I0 > mean2(I0));

centroid = calculateCentroid(I,4);
if(centroid(3) >0.5), weights = flipdim(weights,3); end
weightedI = weights.*I; weightedI = weightedI/max(weightedI(:));
if(centroid(3) >0.5), weightedI = flipdim(weightedI,3); end 

% Find local maximum of surface
[~, maxI] = max(weightedI,[],3);
maxI = padarray(maxI,blockSize*[1 1]/2,'symmetric','both');


maxI_ = im2col(maxI,blockSize*[1 1],'distinct');
maxI_ = max(maxI_);
maxI = reshape(maxI_,sqrt(numel(maxI_))*[1,1]);

[sB{1}, sB{2}] = meshgrid((1:blockSize:blockSize*(size(maxI,1))),(1:blockSize:blockSize*(size(maxI,1))));
% "Surface Fitting using gridfit"
% by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
s4{3} = gridfit(sB{1},sB{2},maxI,s4Idx{1},s4Idx{2},'smoothness',smoothness,...
    'regularizer' ,'laplacian');

if(centroid(3) > 0.5), s4{3} = abs(s4{3}-size(s1Idx{3},2)); end 
s = s4;

end
%%
function C = calculateCentroid(I0,resampleFactor)
sizeI = size(I0);

idx = cell(1,3); dI = cell(1,3);
for i = 1:3, idx{i} = 1:resampleFactor:sizeI(i)-1; end
I = I0(idx{:});
for i = 1:3, idx{i} = idx{i}/sizeI(i); end

[m{1}, m{2}, m{3}] = meshgrid(idx{:});
for i = 1:3, dI{i} = I.*m{i}; end

C = cellfun(@(x) sum(x(:)), dI)/sum(I(:));
end

function [n, tri] = funTriInterpNormals(s,centroid)
% tr:  delauney tiranguluzation of the surface. see TriRep()
% s{1}, s{2}: meshgrid that the tractions are defined on
% fnX, fnY, fnZ: interpolated surface
% inCenter: interpolated center positions of the triRep

dt = DelaunayTri(s{1}(:),s{2}(:));
tri = TriRep(dt.Triangulation,[s{1}(:) s{2}(:) s{3}(:)]);

n_ = faceNormals(tri);
centers = incenters(tri);

nFun = cell(1,3); n = cell(1,3);
for i = 1:3
    nFun{i} = TriScatteredInterp(centers(:,1),centers(:,2),n_(:,i));
    n{i} = nFun{i}(s{1},s{2});
    n{i} = inpaint_nans(n{i});
    if centroid(3) > 0.5, n{i} = -n{i}; end
end

end

function angle = calculateNormalAngle(fn,n)

for i = 1:3, normalVector{i} = n(i)*ones(size(fn{i})); end
fnMag = cellfun(@(x) (x.^2), fn,'UniformOutput',false);
fnMag = sqrt(fnMag{1} + fnMag{2} + fnMag{3});

normMag = cellfun(@(x) (x.^2), normalVector,'UniformOutput',false);
normMag = normMag{1} + normMag{2} + normMag{3};

fnTimesNorm = cellfun(@times, fn, normalVector,'UniformOutput',false);
fnTimesNorm = fnTimesNorm{1} + fnTimesNorm{2} + fnTimesNorm{3};

cosAngle = fnTimesNorm./(normMag.*fnMag);

angle = acosd(cosAngle);

angleMean = mean(angle(:));
angleStd = std(angle(:));
end

function plotSurface(tri,s,n,angle001)
fullScreenFigure();
ha = tight_subplot(2,2);
axes(ha(1));

% surface plot with normals
trisurf(tri); shading flat; hold on;
quiver3(s{1},s{2},s{3},n{1},n{2},n{3});
view(3); axis equal; colorbar; hold off;
title('surface with normal vectors');

axes(ha(2));
I = s{3};
% I = padarray(I,[32,32],'symmetric');
FFT = fft2(I);
FFT = real(FFT);
FFT(1) = nan;
FFT = fftshift(FFT);
imagesc(FFT); axis image; colorbar;
title('FFT of the 3rd surface component');
% FFT = unpadarray(FFT,[32,32],[],'both');

axes(ha(3));
imagesc(angle001); axis image; colorbar;
title('Plot of the 001 angle');

axes(ha(4));
hist(angle001(:),1000);
title('Histogram of the 001 angle');

end

function varargout = parseInputs(varargin)

fileDir = varargin{1};
dataChannel = varargin{2};
dm = varargin{3};
smoothness = varargin{4};
if length(varargin)<5, blockSize = 4;
else blockSize = varargin{5};
end

if length(varargin)<6
    plot01 = 0;
else
    plot01 = varargin{6};
end
% create directory of filenames

if isnumeric(fileDir)
    I = {fileDir};
else
    if ~isstruct(fileDir)
        ext = fileDir(end-3:end);
        if ~strcmpi(ext,'.mat'), fileDir = dir([fileDir,'.mat']);
        else fileDir = dir(fileDir);
        end
    end
    if isempty(fileDir), error('myApp:argChk', 'File name doesn''t exist'); end
    fileNames = {fileDir.name}';
    
    % Load files into memory
    maxTime = length(fileNames);
    I = cell(1,maxTime);
    for i = 1:maxTime,
        I{i} = loadFast(fileNames{i},dataChannel);
    end
end
I = cellfun(@double,I,'UniformOutput',0);

% Apply weights to image
sizeI = size(I{1});
linWeights = shiftdim(linspace(0,1,sizeI(3)),-1);
weights = repmat(linWeights,[sizeI(1),sizeI(2),1]);

if smoothness <= 0, smoothness = 0.0001; end

varargout{1} = I;
varargout{2} = dm;
varargout{3} = weights;
varargout{4} = smoothness;
varargout{5} = blockSize;
varargout{6} = plot01;
end

function vol = loadFast(fileName,channel)
vol = load(fileName);
fName = fieldnames(vol);
vol = getfield(vol,fName{1});
vol = vol{channel};
end