% Plot the background climatological sections from the Canek transects
clear all; close all; clc
%
% Grid resolution (verticla [m] and horizontal [degrees])
vres = 20; hres = 0.03;
% interpolation parameters for both sections
params_climy = [75000, 900, 0.01];
params_climf = [90000, 900, 0.01];

% put your path to the main folder
yourPath = 'C:/Users/nucle/Tesis/Papers/Paper_I/Figures/Programas/CreaRepos/';

% add the external fncs folder to the path:
addpath(genpath([yourPath, '/extrn']));

% path to the files
origin = [yourPath, '/Data/BackgroundSections/'];
% ===============================================================

close all
figure('Position', [70 50 560 604.8889]);
colormap(jet(40))


% Mapping the climatology from the Yucatan section
% load the cell-averages from the Yucatan section
load([origin, '\yuc_vel_clim_clean.mat']);

% make the scatter plot with the normal-to-the-section velocity as color
subplot(2, 2, 1)
box on;
scatter( Nxp(:),Nzp(:), 15,Nvp(:), 'filled' );
hold on;
caxis([-0.1, 1.2]);
set(gca, 'XtickLabels', {})
ylabel('Depth (m)', 'interpreter','latex')
set(gca, 'FontSize', 11);
batsec('yuc', gca, 'k', 2);
axis([-86.8 -84.87, -2100, 100]);

% prepare the Yucatan grid Xy, Zy
[~, ~, ~, xmapy, zmapy, celday, Xy, Zy, dayu]= prep_sec('yuc', vres, hres);

An = Nvp(:) - mean(Nvp(:)); % Remove the overal spatial mean
% Perform the objetctive mapping
[ClimVy, ~, Err]=objan_campomedio_yuc2(An(:),Nxp(:),Nzp(:),params_climy, ...
    xmapy,zmapy);

ClimVy = ClimVy + mean(Nvp(:)); % add-up the spatial mean

% plot the mapped climatological section: Yucatan
Uy = nan(size(Xy));
Uy(celday) = ClimVy(:);


subplot(2, 2, 3)
box on;
pcolor(Xy, Zy, Uy); shading interp
hold on;
[cony, hhy]= contour(Xy, Zy, Uy, [1, 0.9, 0.7, 0.5, 0.3, 0.1, -0.05, 0, -0.1, -0.2], 'k');
clabel(cony, hhy, 'fontsize', 8)
caxis([-0.1, 1.2]);
ylabel('Depth (m)', 'interpreter','latex')
set(gca, 'FontSize', 11);
batsec('yuc', gca, 'k', 2);
axis([-86.8 -84.87, -2100, 100]);
xlabel('Longitude', 'interpreter','latex')

%==============================================================

% Mapping the climatology from the Florida section
% load the cell-averages from the Yucatan section
load([origin, '\flo_vel_clim_clean.mat']);

% make the scatter plot with the normal-to-the-section velocity as color
subplot(2, 2, 2)
box on;
scatter( Nxp(:),Nzp(:), 15,Nvp(:), 'filled' );
hold on;
set(gca, 'XtickLabels', {}, 'YtickLabels', {})
batsec('flo', gca, 'k', 2);
axis([23.1, 24.67, -2100, 100]);
caxis([-0.1, 1.2]);
set(gca, 'FontSize', 11);


% prepare the Florida grid
[mask, xx, yy, xmapf, zmapf, celdaf, Xf, Zf, dafl]= prep_sec('flo', vres, hres);

An = Nvp(:) - mean(Nvp(:)); % Remove the overal spatial mean

[ClimVf, ~, ~]=objan_campomedio_flo2(An(:),Nxp(:),Nzp(:),params_climf,...
    xmapf,zmapf);

ClimVf = ClimVf + mean(Nvp(:)); % add-up the spatial mean


% plot the mapped climatological section: Yucatan

Uf = nan(size(Xf));
Uf(celdaf) = ClimVf(:);

subplot(2, 2, 4)
pcolor(Xf, Zf, Uf);  shading interp
hold on;
[cony, hhy]= contour(Xf, Zf, Uf, [1, 0.9, 0.7, 0.5, 0.3, 0.1, -0.05, 0, -0.1, -0.2], 'k');
clabel(cony, hhy, 'fontsize', 8)
batsec('flo', gca, 'k', 2);
axis([23.1, 24.67, -2100, 100]);
caxis([-0.1, 1.2]);
set(gca, 'YtickLabels', {})


cb = colorbar;
cb.Position = [0.8610 0.1213 0.0311 0.2145];
cb.AxisLocation = 'in';
cb.FontSize = 11;
set(gca, 'FontSize', 11);
ylabel(cb, '(m s$^{-1}$)', 'interpreter','latex', 'fontsize', 13)
xlabel('Latitude', 'interpreter','latex')
%==============================================================

%
set (findobj (gcf, '-property', 'interpreter'), 'interpreter','latex')
set (findobj (gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter','latex');
