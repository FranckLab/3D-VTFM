%% Example run file for the IDVC - LDTFM
% The two images that are going to be run are 'vol00.mat' and 'vol01.mat.'
% the deformation are defined as a four-pole Gaussian prescribed
% displacement field.  See below for distriub
%
% Central z (x_3) plane. Number denots width of gaussian in voxels.  See
% section 3.1 and 3.2 in Bar-Kochba et al. (2014)
% --------------------------
% |                        |
% |     32           64    |
% |                        |
% |                        |
% |                        |
% |     96           128   |
% |                        |
% |                        |
% --------------------------
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   filename: string for the filename prefix for the volumetric images in
%             the current directory.
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'filename*.mat' or 'filename*'
%
%             --- If image is within a cell that contains multichannels ---
%             2) filename{1} = 'filename*.mat' or 'filename*' and
%                filename{2} = channel number containing images you want to
%                              run IDVC on.
%                (if the channel is not provided, i.e. length(filename) = 1
%                , then channel = 1
%
%   sSize: interrogation window (subset) size for the first iterations.
%          Must be 32,64,96, or 128 voxels and a three column
%          array (one for each dimenision) or scalar (equal for all
%          dimensions).
%   incORcum: string that defines the method of running IDVC. Options:
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDVC. Format: cell array,
%      which is a 3D vector (components in x,y,z)  per each time point
%      (units are in voxels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement in z-direction at t=time of size MxNxP
%   cc: peak values of the cross-correlation for each interrogation
%
% NOTES
% -------------------------------------------------------------------------
% To run you need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

clear; close all;

sSize = [96 128 64];
incORcum = 'incremental';
filename = 'vol*.mat';

% Estimate displacements via IDVC
[u, cc, dm] = funIDVC(filename, sSize, incORcum);

%% LDTFM
% fun3DVTFM calculates 3D Traction Force Microscopy metrics
% [Fij, Sij, Eij, ti, tiPN] = fun3DVTFM(u,dm, surface, normals,
% materialProps, materialModel) calculates Deformation gradient Fij,
% Cauchy stress Sij, Lagrangian strain Eij, Strain energy density Uhat,
% surface traction vector ti, and in-plane and out-of-plane
% components of surface tractions
%
% If used please cite:
% Toyjanova J., Hannen, E., Bar-Kochba E., Darling, E.M., Henann, D.L., 
% and Franck, C., (2014) 3D Viscoelastic Traction Force Microscopy.
% Soft Matter 
% doi: 10.1039/c4sm01271b

% Give spacing used in the meshgrid of the results
sizeI = (size(u{1}{1})-1)*dm;

[surface{1}{1},surface{1}{2}] = meshgrid(1:dm:sizeI(1)+1,1:dm:sizeI(2)+1);
surface{1}{3} = 192/2*ones(size(surface{1}{1}))+1;

normals{1}{1} = zeros(size(surface{1}{1}));
normals{1}{2} = zeros(size(surface{1}{1}));
normals{1}{3} = ones(size(surface{1}{1}));

materialModel = 'neohookean'; % material model for the left branck spring 
materialProps = [3600, 0.3, 700, 800];  % [Young's Modulus, Poisson's Ratio, 
                               %  Non-equlibrium Shear modulus (right branch),
                               %  Maxwell Element Viscosity (right branch)]
dt = 100; %seconds
[surface, normals] = calculateSurfaceUi(surface(1), normals(1), u);
% % If you want to calculate the surface based on your pattern (or beads) and
%  the displacements use the following code.
%  Input the name of the volumetric image
%    filePrefix = 'vol*';
%    F = dir(filePrefix);
%    I0 = importdata(F(1).name);
%    dataChannel = 1;
%     [s1,n1,angle001] = findSurface(I{dataChannel}, dataChannel, dm, 12,16,0);
%     [surface, normals] = calculateSurfaceUi(s1, n1, u);
% % If you want to calculate surface normals from already pre-defined
% %  surface grid, you can use the following code:
%  [normals, tri] = calculateNormals(surface)

% Calculate VTFM metrics
[Fij, Sij, Eij, ti, tiPN] = fun3DVTFM(u, dm, surface, normals, materialModel, materialProps, dt);

for i = 1:9
    gradU = sqrt((Fij{1}{1,1}-1).^2 + (Fij{1}{1,2}).^2 + (Fij{1}{1,3}).^2 + ...
        (Fij{1}{2,1}).^2 + (Fij{1}{2,2}-1).^2 + (Fij{1}{2,3}).^2 + ...
        (Fij{1}{3,1}).^2 + (Fij{1}{3,2}).^2 + (Fij{1}{3,3}-1).^2);
end

%% PLOTTING
plotIdx = cell(1,3);
for i = 1:3, plotIdx{i} = 1:dm:sizeI(i)+1; end

close all;
% plot IDVC results
figure;
for i = 1:3
    subplot(2,3,i);
    [~,h] = contourf(plotIdx{1}, plotIdx{2}, u{1}{i}(:,:,ceil(size(u{1}{i},3)/2)),25); colorbar
    set(h,'linestyle','none'); axis image
    title(['displacement component u_',num2str(i)])
    xlabel('X_1'); ylabel('X_2');
    
    subplot(2,3,i+3);
    [~,h] = contourf(plotIdx{1}, plotIdx{3}, squeeze(u{1}{i}(:,48,:))',25); colorbar
    set(h,'linestyle','none'); axis image
    xlabel('X_2'); ylabel('X_3');
    
end

% plot gradU results
figure;
subplot(2,1,1);
[~,h] = contourf(plotIdx{1}, plotIdx{2}, gradU(:,:,ceil(size(gradU,3)/2)),25); colorbar
set(h,'linestyle','none'); axis image
title('|\nablau|')
xlabel('X_1'); ylabel('X_2');

subplot(2,1,2);
[~,h] = contourf(plotIdx{1}, plotIdx{3}, squeeze(gradU(:,48,:))',25); colorbar
set(h,'linestyle','none'); axis image
xlabel('X_2'); ylabel('X_3');

figure;
for i = 1:3
    subplot(1,3,i);
    [~,h] = contourf(plotIdx{1}, plotIdx{2}, ti{1}{i},100); colorbar
    set(h,'linestyle','none'); axis image
    title(['traction component t_',num2str(i)])
    xlabel('X_1'); ylabel('X_2');
end

