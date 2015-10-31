function [Fij, Sij, Eij, Uhat, ti, tiPN] = fun3DTFM(varargin)
% FUNCTIONRUNTFM calculates 3D Traction Force Microscopy metrics
% [Fij, Sij, Eij, Uhat, ti, tiPN] = functionRunTFM(u, dm, surface, normals,
% materialProps, materialModel) calculates Deformation gradient Fij,
% Cauchy stress Sij, Lagrangian strain Eij, Strain energy density Uhat, 
% surface traction vector ti, and in-plane and out-of-plane
% components of surface tractions
%
% INPUTS
% ------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDVC. Format: cell array,
%      which is a 3D vector (components in x,y,z)  per each time point
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement in z-direction at t=time of size MxNxP
%   dm: the spacing of the meshgrid generated in FIDVC
%   surface: cell array containing the surface coordinate grid defined by a
%            rectangular meshgrid (can be non-uniform)
%            surface{time}{1} = x-coordinates of nodes
%            surface{time}{2} = y-coordinates of nodes
%            surface{time}{3} = z-coordinates of nodes
%   normals: cell array containing the normal vectors for each point on the 
%            surface
%            normals{time}{1} = x-component of normal vector
%            normals{time}{2} = y-component of normal vector
%            normals{time}{3} = z-component of normal vector
%   materialModel: 'neohookean', 'linearElastic'; Default = 'linearElastic'
%   materialProps: material Properties of the substrate.
%                  Array: [Young's Modulus, Poisson's ratio]
% 
% OUTPUTS
% ------------------------------------------------------------------------
%   For all tensor outputs, Fij, Sij, and Eij we adapt the following
%   notation. Consider a second-order tensor Aij with nine-components for
%   each time increments. Each component is a cell located within another
%   cell that denotes time,
%                    |{A11}, {A12}, {A13}|
%        Aij{time} = |{A21}, {A22}, {A23}|
%                    |{A11}, {A12}, {A13}|
%   For example, to access the A12 for the 1st time increment you would
%   call Aij{1}{1,2}.
%   
%   Variable definitions:
%       Fij: deformation gradient tensor
%       Sij: Cauchy stress tensor
%       Eij: Lagrangian Strain Tensor (neohookean) or Infinitesimal strain
%            tensor (linear elastic)
%       Uhat: strain energy density
%       ti: tractions in the global coordinate system
%             ti{time}{1} = x-component
%             ti{time}{2} = y-component
%             ti{time}{3} = z-component
%             ti{time}{4} = magnitude
%       tiPN: tractions in the local coordinate system
%             tiPN{time}{1}: parallel, shear, or in-plane component
%             tiPN{time}{2}: perpendicular, normal, or out-of-plane component
%        
% NOTES
% ------------------------------------------------------------------------
% The function was tested on 2013b Matlab. Please note that certain
% triangulation functions on older or newer versions may not be
% compatible mostly due to syntax.
%
% 
% If used please cite:
% Toyjanova J., Bar-Kochba E., López-Fagundo C., Reichner J., 
% Hoffman-Kim D, Franck, C. (2014) High Resolution, Large Deformation 3D
% Traction Force Microscopy. PLoS ONE 9(4): e90976. 
% doi:10.1371/journal.pone.0090976

[u, dm, surface, normals, materialModel, materialProps, OS] ...
    = parseInputsTFM(varargin{:})

% DEFORMATION GRADIENT
[Fij, J] = calculateFij(u,dm,'optimal9');

% STRESS TENSOR
[Sij, Eij, Uhat] = calculateSij(Fij,J,materialProps,materialModel);

% CALCULATE SURFACE TRACTIONS (Cauchy Relation)
[ti, tiPN] = calculateTractions(Sij,surface,normals,OS);

end

%% ========================================================================
function  varargout = parseInputsTFM(varargin)
u = varargin{1};
dm = varargin{2};
surface = varargin{3};
normals = varargin{4};
materialProps = varargin{6};

if isempty(varargin{5})
    materialModel = 'linearElastic';
else
    materialModel = varargin{5};
end

if strcmpi(materialModel,'linearElastic'), OS = 0; % reference configuration
else OS = 1; % deformed configuration
end

varargout{      1} = u;
varargout{end + 1} = dm;
varargout{end + 1} = surface;
varargout{end + 1} = normals;
varargout{end + 1} = materialModel;
varargout{end + 1} = materialProps;
varargout{end + 1} = OS;

end

%% ========================================================================
function [Sij, Eij, U] = calculateSij(Fij, J, materialProps, materialModel)
% Stress Based on Isotropic, linear elastic material behavior
% [Sij,Eij,Uhat] = calculateSij(Fij, J, [E, nu], type)
%
% INPUT: Fij, J, materialProps, type
% -------------------------------------------------------------------------
% Fij: deformation gradient tensor
% J: Jacobian matrix
% materialProps: vector [E,nu], where E: Young's Modulus and nu: Poisson's
%                ratio
% type: 'neohookean', 'linearElastic'
%
% OUTPUT: Sij, Eij, Uhat
% -------------------------------------------------------------------------
% Sij: Cauchy stress tensor
% Eij: strain tensor (infinitesimal or lagrangian)
% Uhat: strain energy density

E = materialProps(1); %Pascals
nu = materialProps(2);

maxTime = length(J);
Sij = cell(maxTime,1);
Eij = cell(maxTime,1);
U = cell(maxTime,1);

for i = 1:maxTime
    switch lower(materialModel)
        case 'neohookean', [Sij{i}, Eij{i}, U{i}] = neoHookean(Fij{i},J{i},E,nu);
        otherwise,         [Sij{i}, Eij{i}, U{i}] = linearElastic(Fij{i},E,nu);
    end
end

end

%% ========================================================================
function [Sij, Eij, U] = linearElastic(Fij,E,nu)
% Linear Elastic Isotropic Solid small deformation
% Bower AF (2011) Applied Mechanics of Solids. CRC press., pg 71

Eij = funCalculateeij(Fij);
Sij = cell(3,3);
trace = Eij{1,1} +  Eij{2,2} +  Eij{3,3};

for i = 1:3
    for j = 1:3
        dij = i==j;
        Sij{i,j} = E/(1+nu)*(Eij{i,j} + nu/(1-2*nu) * trace * dij);
    end
end

% Strain energy density
U = 0;
for k = 1:9
    U = U + Sij{k}.*Eij{k};
end
U = U/2;

end

%% ========================================================================
function [Sij, Eij, U] = neoHookean(Fij,J,E,nu)
% Generalized neo-Hookean solid (adapted from Treloar 1948)
% Bower AF (2011) Applied Mechanics of Solids. CRC press., pg 100

Bij = funCalculateBij(Fij);
Sij = cell(3,3);

mu = E/(2*(1+nu));
K = E/(3*(1-2*nu));

muJ = mu/J.^(5/3);
KJ = K*(J-1);

trB = 1/3*(Bij{1,1} + Bij{2,2}+ Bij{3,3});
for i = 1:3
    for j = 1:3
        if i == j, Sij{i,j} = muJ.*(Bij{i,j} - trB) + KJ;
        else Sij{i,j} = muJ.*Bij{i,j}; end
    end
end

% Strain energy density
I1bar = 3*trB./J.^(2/3);
U = mu/2*(I1bar - 3) + K/2*(J - 1).^2;

Eij = funCalculateLagrangianEij(Fij);
end

%% ========================================================================
function Eij = funCalculateeij(Fij)
% Bower AF (2011) Applied Mechanics of Solids. CRC press., pg 22
Eij = cell(3,3);
Eij{1,1} = Fij{1,1} - 1;
Eij{2,2} = Fij{2,2} - 1;
Eij{3,3} = Fij{3,3} - 1;

Eij{1,2} = 0.5*(Fij{1,2} + Fij{2,1});
Eij{1,3} = 0.5*(Fij{1,3} + Fij{3,1});
Eij{2,3} = 0.5*(Fij{2,3} + Fij{3,2});

Eij{2,1} = Eij{1,2}; Eij{3,1} = Eij{1,3}; Eij{3,2} = Eij{2,3};

end

%% ========================================================================
function Eij = funCalculateLagrangianEij(Fij)
Cij = num2cell(zeros(3));

for i = 1:3
    for j = 1:3
        for k = 1:3, Cij{i,j} = Cij{i,j} + Fij{k,i}.*Fij{k,j}; end
    end
end
% Bower AF (2011) Applied Mechanics of Solids. CRC press., pg 20
Eij = cell(3,3);
for i = 1:3
    for j = 1:3
        I = double(i==j);
        Eij{i,j} = 0.5*(Cij{i,j} - I);
    end
end

end

%% ========================================================================
function  Bij = funCalculateBij(Fij)
% Bij = funCalculateBij(Fij)
% Calculates Cauchy-Green Deformation Tensor Bij = Fik*Fjk
% Bower AF (2011) Applied Mechanics of Solids. CRC press., pg 96
Bij = num2cell(zeros(3));

for i = 1:3
    for j = 1:3
        for k = 1:3, Bij{i,j} = Bij{i,j} + Fij{i,k}.*Fij{j,k}; end
    end
end

end

%% ========================================================================
function [Fij, J, gradU] = calculateFij(varargin)
% [Fij, J] = calculateFij(u,spacing,type {'optimal9'})
% [Fij, J, gradU] = calculateFij(u, dm, type)
%
% INPUT: u, dm, type
% ------------------------------------------------------------------------
% u: displacement vector from FIDVC
% dm: spacing of the meshgrid
% type: is a string with a kernel choice.
%       Options:'prewitt', 'sobel','scharr', 'stencil', 'optimal5' - 'optimal19';
%
% OUTPUT: Fij, J, gradU
% ------------------------------------------------------------------------
% Fij: deformation gradient tensor
% J: Jacobian matrix
% gradU: gradient of displacement vector.
% The format of each of the outputs is cell arrays with
% Data{time}{component}. Note, Fij is a deformation gradient tensor with 9
% components.
%
% NOTES
% ------------------------------------------------------------------------
% For more information on optimal tap filter differentiation please
% see:
% Farid H, Simoncelli EP (2004) "Differentiation of discrete multidimensional
% signals." IEEE transactions on image processing : a publication of the IEEE
% Signal Processing Society 13: 496–508. doi: 10.1109/tip.2004.823819

[u,spacing,type] = parseInputs(varargin{:});

maxTime = length(u);
Fij = cell(maxTime,1);
J = cell(maxTime,1);
gradU = cell(maxTime,1);

for i = 1:maxTime
    Fij{i} = funCalculateFij(u{i},spacing,type);
    J{i} = funCalculateJ(Fij{i});
    gradU{i} = funCalculateGradU(Fij{i});
    
end

end

%% ========================================================================
function Fij = funCalculateFij(u,spacing,type)
% FUNCALCULATEFIJ calculates Displacement Gradient
Fij = cell(3,3);
for i = 1:3,
    [Fij{i,1}, Fij{i,2}, Fij{i,3}] = gradientN(u{i},type);
end

for k = 1:9, Fij{k} = Fij{k}/spacing; end
for i = 1:3, Fij{i,i} = Fij{i,i} + 1; end

end

%% ========================================================================
function J = funCalculateJ(Fij)
% FUNCALCULTATEJ  calculates Jacobian matrix J = det(F)
J =     Fij{1,1}.*Fij{2,2}.*Fij{3,3};
J = J + Fij{1,2}.*Fij{2,3}.*Fij{3,1};
J = J + Fij{1,3}.*Fij{2,1}.*Fij{3,2};
J = J - Fij{1,3}.*Fij{2,2}.*Fij{3,1};
J = J - Fij{1,2}.*Fij{2,1}.*Fij{3,3};
J = J - Fij{1,1}.*Fij{2,3}.*Fij{3,2};
end

%% ========================================================================
function gradU = funCalculateGradU(Fij)

for i = 1:9
    gradU = sqrt((Fij{1,1}-1).^2 + (Fij{1,2}).^2 + (Fij{1,3}).^2 + ...
        (Fij{2,1}).^2 + (Fij{2,2}-1).^2 + (Fij{2,3}).^2 + ...
        (Fij{3,1}).^2 + (Fij{3,2}).^2 + (Fij{3,3}-1).^2);
end

end

%% ========================================================================
function [u,spacing,type] = parseInputs(varargin)
u = varargin{1};
spacing = varargin{2};

if length(varargin) < 3, type = 'optimal9';
else type = varargin{3};
end

for i = 1:length(u)
    u{i} = cellfun(@double, u{i}, 'UniformOutput', 0);
end

end

%% ========================================================================
function varargout = gradientN(varargin)
%    [Vx, Vy, Vz] = gradientN(V,type)
[I, type] = parseInputsGN(varargin{:});
[w,d] = fspecialSeperable(type);
padOpt = 'symmetric';

if ismatrix(I)
    if length(w) > 1
        Ix = -imfilter(I,w','same',padOpt);
        Iy = -imfilter(I,w,'same',padOpt);
    else
        Ix = -I; Iy = -I;
    end
    
    Ix = -imfilter(Ix,d,'conv','same',padOpt);
    Iy = -imfilter(Iy,d','conv','same',padOpt);
else
    if length(w) > 1
        Ix = -imfilter(I,w','same',padOpt); Ix = -imfilter(Ix,reshape(w,1,1,[]),'same',padOpt);
        Iy = -imfilter(I,w,'same',padOpt);  Iy = -imfilter(Iy,reshape(w,1,1,[]),'same',padOpt);
        Iz = -imfilter(I,w,'same',padOpt);  Iz = -imfilter(Iz,w','same',padOpt);
    else
        Ix = -I; Iy = -I; Iz = -I;
    end
    
    Ix = -imfilter(Ix,d,'same',padOpt);
    Iy = -imfilter(Iy,d','same',padOpt);
    Iz = -imfilter(Iz,reshape(d,1,1,[]),'same',padOpt);
    varargout{3} = Iz;
end

varargout{1} = Ix; varargout{2} = Iy;
end

%% ========================================================================
function [w,d] = fspecialSeperable(type)
sumFlag01 = 1;
switch lower(type)
    case 'fb', w = 1; d = [1 -1];
    case 'prewitt', w = [1 1 1]; d = [1 0 -1];
    case 'sobel', w = [1 2 1]; d = [1 0 -1];
    case 'scharr', w = [3 10 3]; d = [1 0 -1];
    case 'stencil', w = 1; d = [-1 8 0 -8 1]/12; sumFlag01 = 0;
    otherwise
        if strcmpi(type(1:7),'optimal')
            sumFlag01 = 0;
            switch str2double(type(8:end))
                case 5,
                    w = [0.0376593171958139,0.249153396177344,0.426374573253684,0.249153396177344,0.0376593171958139];
                    d = [0.109603762960256,0.276690988455550,0,-0.276690988455550,-0.109603762960256];
                case 7,
                    w = [0.00541196740974425,0.0695905825057286,0.244559723794791,0.360875452579473,0.244559723794791,0.0695905825057286,0.00541196740974425];
                    d = [0.0194786630434688,0.123914729925000,0.193554838845068,0,-0.193554838845068,-0.123914729925000,-0.0194786630434688];
                case 9,
                    w = [0.000737598362084457,0.0155298478667911,0.0902598182960924,0.234469365350285,0.318006740249494,0.234469365350285,0.0902598182960924,0.0155298478667911,0.000737598362084457];
                    d = [0.00303163095459785,0.0352414678518254,0.118879484725614,0.144382960330377,0,-0.144382960330377,-0.118879484725614,-0.0352414678518254,-0.00303163095459785];
                case 11,
                    w = [9.73660046309032e-05,0.00304151665483548,0.0261779198435670,0.103249150042758,0.223724696218449,0.287418702471518,0.223724696218449,0.103249150042758,0.0261779198435670,0.00304151665483548,9.73660046309032e-05];
                    d = [0.000439905393637534,0.00808328124581577,0.0449782033510644,0.108830913335368,0.112870645898035,0,-0.112870645898035,-0.108830913335368,-0.0449782033510644,-0.00808328124581577,-0.000439905393637534];
                case 13,
                    w = [1.49541584976590e-05,0.000608044910058278,0.00691160005872185,0.0369033951658625,0.112173742521741,0.212488560211325,0.261799405947588,0.212488560211325,0.112173742521741,0.0369033951658625,0.00691160005872185,0.000608044910058278,1.49541584976590e-05];
                    d = [7.14573394720322e-05,0.00177958343198347,0.0138717416627181,0.0506584387298514,0.0970421434222207,0.0891267488305157,0,-0.0891267488305157,-0.0970421434222207,-0.0506584387298514,-0.0138717416627181,-0.00177958343198347,-7.14573394720322e-05];
                case 15,
                    w = [5.98860795150199e-06,0.000191445809195946,0.00198513495426173,0.0115248319580300,0.0440836748027682,0.114823905477918,0.203925200191027,0.246919636397694,0.203925200191027,0.114823905477918,0.0440836748027682,0.0115248319580300,0.00198513495426173,0.000191445809195946,5.98860795150199e-06];
                    d = [2.68289243388421e-05,0.000526939396634328,0.00397615629561318,0.0177322281254572,0.0506284057542549,0.0879544206261874,0.0780467781003990,0,-0.0780467781003990,-0.0879544206261874,-0.0506284057542549,-0.0177322281254572,-0.00397615629561318,-0.000526939396634328,-2.68289243388421e-05];
                case 17,
                    w = [7.69392199800271e-07,3.43462804426615e-05,0.000460563932733138,0.00328564897335160,0.0153609347354413,0.0504415385506543,0.117907645586470,0.196238466137000,0.232540172823407,0.196238466137000,0.117907645586470,0.0504415385506543,0.0153609347354413,0.00328564897335160,0.000460563932733138,3.43462804426615e-05,7.69392199800271e-07];
                    d = [3.73783393466673e-06,0.000104597678828706,0.00102618747708425,0.00569710286194269,0.0209031813657965,0.0513851162956773,0.0801001019942813,0.0666269955155201,0,-0.0666269955155201,-0.0801001019942813,-0.0513851162956773,-0.0209031813657965,-0.00569710286194269,-0.00102618747708425,-0.000104597678828706,-3.73783393466673e-06];
                case 19,
                    w = [1.52643475621529e-07,7.90799631627647e-06,0.000122710209565798,0.000998814542578885,0.00531312830085921,0.0202775623584220,0.0572293928770952,0.120087013816114,0.187334615968446,0.217257402574254,0.187334615968446,0.120087013816114,0.0572293928770952,0.0202775623584220,0.00531312830085921,0.000998814542578885,0.000122710209565798,7.90799631627647e-06,1.52643475621529e-07];
                    d = [7.63266760902856e-07,2.52762743204803e-05,0.000290470367850947,0.00185854511355811,0.00795173988567152,0.0240580568530862,0.0508916097992498,0.0712071454435055,0.0555263098864427,0,-0.0555263098864427,-0.0712071454435055,-0.0508916097992498,-0.0240580568530862,-0.00795173988567152,-0.00185854511355811,-0.000290470367850947,-2.52762743204803e-05,-7.63266760902856e-07];
                otherwise
                    w = [0.229878817299031,0.540242365401939,0.229878817299031];
                    d = [0.425286806086887,0,-0.425286806086887];
            end
        else
            w = 1; d = [-1 0 1];
        end
end

w = w/sum(abs(w));

if sumFlag01
    d = d/sum(abs(d));
end


end

%% ========================================================================
function [I, type] = parseInputsGN(varargin)

if length(varargin) < 2, varargin{2} = 'prewitt'; end

I = double(varargin{1});
type = varargin{2};
end

%% ========================================================================
function [ti, tiPN] = calculateTractions(varargin)
% [ti, tiPN] = calculateTractions(Sij,s,n,OS); is used to calculate
% surface traction 3D vector
%
% INPUT: Sij, s, n, OS
% ------------------------------------------------------------------------
% Sij: Cauchy Stress tensor
% s: surface coordinate grid
% n: surface normals
% OS: 1 -  current configuration (Hyperelastic model) or 0 - reference
%           configuration (Linear elastic model)
%
% OUTPUT: ti, tiPN
% -----------------------------------------------------------------------
% ti: surface traction 3D vector components per each time stack.
%     Output format - ti{timePoint}{component}. ti{t}{4} - magnitude
% tiPN: in-plane and out-of-plane components of surface traction vector per
%       each time stack.
%       Output format: tiPN{t}{1}: In-plane or tiPN{t}{2}: Out-of-plane}
%
% NOTES
% ----------------------------------------------------------------------
% For more information please see Section "Estimating 3D LD Cellular
% Tractions." in Toyjanova et al.,2014(DOI:10.1371/journal.pone.0090976)

[Sij,m,s,n,OS] = parseInputsTi(varargin{:});

interpMethod = 'cubic';
maxTime = length(Sij);

ti = cell(maxTime,1);
for i = 1:maxTime
    ti{i} = funCalculateTractions(Sij{i},m,s{i+OS},n{i+OS},interpMethod);
end


tiPN = cell(maxTime,1);
for k = 1:maxTime
    tiNormal{k} = abs(ti{k}{1}.*n{k+1,1}{1}+ti{k}{2}.*n{k+1,1}{2}+ti{k}{3}.*n{k+1,1}{3});
    for m = 1:3
        tiParallel_{k}{m} = ti{k}{m}-tiNormal{k}.*n{k+1,1}{m};
    end
    tiParallel{k,1} = sqrt(tiParallel_{k}{1}.^2+tiParallel_{k}{2}.^2+tiParallel_{k}{3}.^2);
    tiPN{k}{1} = tiParallel{k};
    tiPN{k}{2} = tiNormal{k};
end


end

%% ========================================================================
function ti = funCalculateTractions(Sij,m,s,n,interpMethod)

s33 = sampleData3(Sij{3,3},m,s,interpMethod);
s11 = sampleData3(Sij{1,1},m,s,interpMethod);
s22 = sampleData3(Sij{2,2},m,s,interpMethod);

s13 = sampleData3(Sij{1,3},m,s,interpMethod);
s23 = sampleData3(Sij{2,3},m,s,interpMethod);
s12 = sampleData3(Sij{1,2},m,s,interpMethod);

% Time to calculate those tractions
ti{1} = s11.*n{1} + s12.*n{2} + s13.*n{3};
ti{2} = s12.*n{1} + s22.*n{2} + s23.*n{3};
ti{3} = s13.*n{1} + s23.*n{2} + s33.*n{3};

for i = 1:3
    ti{i} = inpaint_nans(ti{i});
    ti{i} = real(ti{i});
end

ti{4} = sqrt(ti{1}.^2 + ti{2}.^2 + ti{3}.^2);


end

%% ========================================================================
function f1 = sampleData3(f0,m,s,method)
f1 = interp3(m{1},m{2},m{3},f0,s{1},s{2},s{3},method);
end

%% ========================================================================
function varargout = parseInputsTi(varargin)
%  [Sij,m,s,n,OS] = parseInputs(Sij,dm,s,n,OS)
Sij = varargin{1};
s = varargin{2};
n = varargin{3};
OS = varargin{4};

mSize = size(Sij{1}{1});
dm = s{1}{2}(2) - s{1}{2}(1);
idx = cell(1,3);
for i = 1:3; idx{i} = (1:dm:mSize(i)*dm); end
[m{1}, m{2}, m{3}] = meshgrid(idx{:});

varargout{      1} = Sij;
varargout{end + 1} = m;
varargout{end + 1} = s;
varargout{end + 1} = n;
varargout{end + 1} = OS;

end