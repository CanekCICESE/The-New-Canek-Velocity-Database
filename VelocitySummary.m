% Used to create the results on the transport time series analysis show in
% the subsection 3a, Figure 5.
% This is also an example code on how to read the main database Netcdf
% G. Durante, 2024
clear all; close all; clc;

% path to the souce NetCDF
% put your path to the main folder
yourPath = 'C:/Users/nucle/Tesis/Papers/Paper_I/Figures/Programas/CreaRepos/';

% add the external fncs folder to the path:
addpath(genpath([yourPath, '/extrn']));

% path to the souce NetCDF
origin = [yourPath, 'Data/'];

arch = [origin, 'New_Canek_Database_2024.nc'];

tid =ncread(arch, 'time'); % read the time vector
lon =ncread(arch, 'ylon'); % read the time longitude
z =ncread(arch, 'ydepth'); % read the time depth

u =ncread(arch, 'Uy'); % read the time depth
v =ncread(arch, 'Vy'); % read the time depth

% read the batimetric profile
xbat =ncread(arch, 'x_bat'); % read the time depth
bat =ncread(arch, 'bat'); % read the time depth

% The NetCDF file was creating by de-rotating the mapped velocities,since
% the volume-conservation constraint was applyied on the along-channel
% velocity.
% however, it is quite the same as the Figure 5.



% rotating u and v to obtain the normal-to and parallel-to-the section
% componets
vn = u;
vm = nanmean(v, 3);
cel = find(~isnan(vm)); % get the indices of non nan values
cm = nan([length(cel), length(tid)]);

for k = 1 : length(tid)
    udum = u(:, :, k); vdum = v(:, :, k);
    udum = udum(cel); vdum = vdum(cel);
    [cn, cp] = rotaVels('yuc', udum(:), vdum(:));
    cm(:, k) = cn(:);
end
%==========================================
% compute the time mean
Ui = nan(size(v, 1, 2));
Ui(cel) = nanmean(cm, 2);

close all

figure;
subplot(1, 2, 1)
chy = pcolor(lon, z, Ui');shading interp;
hold on;
[cony, hhy]= contour(lon, z, Ui', [1, 0.9, 0.7, 0.5, 0.3, 0.1, -0.05, 0, -0.1, -0.2], 'k');
clabel(cony, hhy)
colormap(jet(40));

plot(xbat, bat, 'color', [1 1 1], 'linewidth', 12);
plot(xbat, bat, 'color', 0.5*[1 1 1], 'linewidth', 2);

ylim([-2100, 100])
axis([-86.8 -84.87, -2100, 100]);
caxis([-0.1, 1.2]);
cb = colorbar;

xlim([-86.8  -84.87])

ylabel('Depth (m)', 'interpreter','latex')

set(gca, 'FontSize', 11);
set (findobj (gcf, '-property', 'interpreter'), 'interpreter','latex')
set (findobj (gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter','latex');

text(-86.74165, -1954, '\textbf{(A)}', 'interpreter','latex', 'fontsize', 14);

tit = title('\textbf{Mean velocity sections (2012-2020)}', 'interpreter','latex');
set(gca, 'TickLength', 0.01*[1 1],...
    'Layer', 'Top', 'Box', 'on', 'YMinorTick',...
    'On', 'XMinorTick', 'On');



%% Make the plots of the Figure 5 (only Yucatan)
subplot(1, 2, 2)

% compute the time std map
Ui(cel) = nanstd(cm, [], 2);
chy = pcolor(lon, z, Ui');shading interp;
hold on;
[cony, hhy]= contour(lon, z, Ui', [0.05, 0.1, 0.15, 0.2, 0.25], 'k');
clabel(cony, hhy)
colormap(jet(40));
clabel(cony, hhy, 'fontsize', 8)


batsec('yuc', gca, rgb('white'), 8);
batsec('yuc', gca, rgb('Gray'), 2);

axis([-86.8 -84.87, -2100, 100]);
caxis([0.001, 0.4]);
ylim([-2100, 100])


set(gca, 'FontSize', 11);
set (findobj (gcf, '-property', 'interpreter'), 'interpreter','latex')
set (findobj (gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter','latex');
xlim([-86.8  -84.87])
text(-86.74165, -1954, '\textbf{(B)}', 'interpreter','latex', 'fontsize', 14);
set(gca, 'TickLength', 0.01*[1 1],...
    'Layer', 'Top', 'Box', 'on', 'YMinorTick',...
    'On', 'XMinorTick', 'On');


colorbar;
title('\textbf{Standard deviation (2012-2020)}', 'interpreter','latex');


