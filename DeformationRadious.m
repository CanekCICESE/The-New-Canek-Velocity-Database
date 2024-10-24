% Estimations of the first baroclinic radious of deformation following (Chelton et al., 1998)
% This code uses the function gsw_Nsquared.m TESO-10 library 
% you can download the library from: https://www.teos-10.org/software/gsw_matlab_v3_06_16.zip

clear all; clc; close all;


% path to the file
yourPath = 'C:/Users/nucle/Tesis/Papers/Paper_I/Figures/Programas/CreaRepos/';

addpath(genpath([yourPath, '\extrn']));

origin = [yourPath, '/Data/CTD/'];

% load the CTD data
load([origin, 'CTD_cnk_Hist.mat'])



Ro = [];
lat = [];
lon = [];
depth = [];

lances = unique(fe); % find unique dates 
tags = ones(size(fe)); % this vector will have the indices or "tags" for each cast
ncasts = length(lances); % number of single casts
aa = 1;
m = 2;

for i = 1 : ncasts
    index = (fe == lances(i)); % find positions of a given cast
    tags(index) = i; % tag the given cast (they are sorted by date) 
    %  Casts are labeled as 1, 2, 3 ... to "ncasts"
    
    sdu = sa(index); 
    tdu = ct(index);
    pdu = p(index);
    ladu = lam(index);
    inba = isnan(sdu + tdu + pdu + ladu);
    
    sdu(inba) = []; tdu(inba) = []; pdu(inba) = []; ladu(inba) = [];
    % determine if the casts are "complete" in depth
    if min(p(index)) <= 40 &  max(p(index)) >= 650  & max(p(index)) < 2100
        
    % interpolate every cast, just in case is not at the same vertical resolution   
    % normally do
    [N2, p_m] = gsw_Nsquared(sdu, tdu, pdu, ladu);
    N2(N2<0) = 0;
    N = sqrt(N2)*2*pi;
    x1 = N; x2 = p_m;
   
    dz = 10;
    x3 = (40:dz:max(p_m))';
    x4 = interp1(x2, x1, x3);
    %=============================================================================
    % compute the Ro using the equation from (Chelton et al. 1998)
    Ro(aa) = ( 1./( coriolis(nanmean(lam(:)))*m*pi ) * nansum(x4)*dz )  * 1e-3;
    
    lat(aa) = mean(lam(index));
    lon(aa) = mean(lom(index));
    depth(aa) = max(p_m);
    mN(aa) = nansum(N);
    dzz(aa) = dz;
    aa = aa + 1;
    
    end
end


%%
close all

% figure;set(gcf, 'pos', [50 50 1414 738])
% axes('pos', [0.5552 0.2290 0.3498 0.6960])
[xv, yv] = muadro([-86.6296 -84.8879], [21.4141 22.0294]);
inyu = inpolygon(lon, lat, xv, yv);

plot(depth(inyu), Ro(inyu), 'ok', 'markerfacecolor',  rgb('RoyalBlue'))

b = polyfit(depth(inyu), Ro(inyu), 3);
yp = polyval(b, 600:2100);
hold on;
plot(600:2100, yp, 'b', 'linewidth', 3)

rmyu = mean(yp)
[min( Ro(inyu)) max( Ro(inyu))]


%

hold on;
[xv, yv] = muadro([-82.1774 -80.4355], [22.9898 24.5000]);
infl = inpolygon(lon, lat, xv, yv);

plot(depth(infl),Ro(infl), 'sk', 'markerfacecolor', rgb('OrangeRed'), 'markersize', 9)

b = polyfit(depth(infl), Ro(infl), 2);
yp = polyval(b, 750:1600);
hold on;
plot(750:1600, yp, 'r', 'linewidth', 3)

rmflo = mean(yp)
[min( Ro(infl)) max( Ro(infl))]

title( '\textbf{Rossby Radius of Deformation (Chelton et al., 1998)}' );

xlabel('Profile depth (m)')
ylabel('Rossby radius of deformation (km)')


ax0(1) = plot(NaN,NaN, 'sk', 'markerfacecolor', rgb('OrangeRed'), 'markersize', 9);
ax0(2) = plot(NaN,NaN, 'ok', 'markerfacecolor',  rgb('RoyalBlue'));


ll = legend(ax0, {'Florida', 'Yucatan'}, 'Location', 'Best');


Fs = findall(gcf,'-property','FontSize');
for k = 1 : length(Fs)
    if strcmp(get(Fs(k), 'type'), 'axes')
        if Fs(k).FontSize < 12
            Fs(k).FontSize = 12;
        end
    end
end

set(findall(gcf,'-property','Interpreter'),'Interpreter', 'latex');
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex');

text(2078, 98.5564, '\textbf{(I)}','Interpreter', 'latex', 'fontsize', 17);

%%
% axes('pos', [0.753889674681754 0.234540527223454 0.147807705924052 0.19308951974927]);
% box on;
% 
% 
% plot(lon(infl), lat(infl), '.', 'color', rgb('OrangeRed'), 'markersize', 7);
% hold on;
% 
% plot(lon(inyu), lat(inyu), '.', 'color',  rgb('RoyalBlue'), 'markersize', 7);
% hold on;
% 
% % pongolfo(gca, rgb('black'), 0);
% axis equal tight
% axis([-87.3, -80.7, 20.5, 25])
% text(-87.1071, 24.2678, '\textbf{(II)}','Interpreter', 'latex', 'fontsize', 17);
% set(gca, 'YtickLabels', {}, 'XtickLabels', {})
