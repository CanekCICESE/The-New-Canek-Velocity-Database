% Used to create the results on the transport time series analysis show in
% the subsection 3a, Volume transport estimations, Table I and Figure 7
%
% G. Durante, 2024

clear all; close all; clc;


% put your path to the main folder
yourPath = 'C:/Users/nucle/Tesis/Papers/Paper_I/Figures/Programas/CreaRepos/';

% add the external fncs folder to the path:
addpath(genpath([yourPath, '/extrn']));

% path to the souce NetCDF
origin = [yourPath, 'Data/'];

arch = [origin, 'New_Canek_Database_2024.nc'];

tid =ncread(arch, 'time'); % read the time vector 

Ty = ncread(arch, 'Ty'); % read the transport in the Yucatan section
Py = ncread(arch, 'Py'); % read the uncertainty of the transport in the Yuc sec.

Tf = ncread(arch, 'Tf'); % read the transport in the Florida section
Pf = ncread(arch, 'Pf'); % read the uncertainty of the transport in the Flo sec.

% remove the NaNs
Tiy = Ty;
diy = tid;
Piy = Py;
Piy(isnan(Tiy)) = [];
diy(isnan(Tiy)) = [];
Tiy(isnan(Tiy)) = [];

Tif = Tf;
dif = tid;
Pif = Pf;
Pif(isnan(Tif)) = [];
dif(isnan(Tif)) = [];
Tif(isnan(Tif)) = [];

% Compute effective degress of fredom and standard error of de mean
[EDF, DE, flgs] = GEL( Tiy, Tiy );
DE = DE(2, 2);

ESy = (round(std(Tiy)/sqrt(DE), 2, 'decimal')); % standard error of the mean Yucatan transp


[EDF, DE, flgs] = GEL( Tif, Tif );
DE = DE(2, 2);
ESf = (round(std(Tif)/sqrt(DE), 2, 'decimal')); % standard error of the mean Florida transp

%%
Quantity = {'Mean (Sv)';'std (Sv)';'Trend slope (Sv/year)';'Annual Amp. (Sv)';'Semi-annual amp. (Sv)';'Bi harm. exp.Var.'};
ColumnNames = {'Quantity', 'Yucatan section', 'Florida section'};

% trend analysis
[estimate, se, ci_margin, weighted_mean] = weighthedRegres( diy(:), Tiy(:), Piy(:) );
% harmnic analysis
[A,Fase,epsA,epsFase,dta,VE,as,epsas]=ajustew(diy,Tiy,2,365.25);


Yuc_Sec = {[num2str(mean(Tiy), 3), ' \pm ', num2str(ESy, 2)]; [num2str(std(Tiy), 2)]; ...
    [num2str(estimate(2)*365.25, 1), ' \pm ', num2str(ci_margin(2)*365.25, 1)]; 
    [num2str(A(2), 1), ' \pm ', num2str(epsA(2), 1)];
    [num2str(A(3), 1), ' \pm ', num2str(epsA(3), 1)];
    [num2str(VE*100, 3), '%']};

% trend analysis
[estimate, se, ci_margin, weighted_mean] = weighthedRegres( diy(:), Tiy(:), Piy(:) );
% harmnic analysis
[A,Fase,epsA,epsFase,dta,VE,as,epsas]=ajustew(diy,Tiy,2,365.25);
% creat coloums for the table Yucatan
Yuc_Sec = {[num2str(mean(Tiy), 3), ' \pm ', num2str(ESy, 2)]; [num2str(std(Tiy), 2)]; ...
    [num2str(estimate(2)*365.25, 1), ' \pm ', num2str(ci_margin(2)*365.25, 1)]; 
    [num2str(A(2), 1), ' \pm ', num2str(epsA(2), 1)];
    [num2str(A(3), 1), ' \pm ', num2str(epsA(3), 1)];
    [num2str(VE*100, 3), '%']};

% trend analysis
[estimate, se, ci_margin, weighted_mean] = weighthedRegres( dif(:), Tif(:), Pif(:) );
% harmnic analysis
[A,Fase,epsA,epsFase,dta,VE,as,epsas]=ajustew(dif,Tif,2,365.25);
% creat coloums for the table Yucatan
Flo_Sec = {[num2str(mean(Tif), 3), ' \pm ', num2str(ESf, 2)]; [num2str(std(Tif), 2)]; ...
    [num2str(estimate(2)*365.25, 1), ' \pm ', num2str(ci_margin(2)*365.25, 1)]; 
    [num2str(A(2), 1), ' \pm ', num2str(epsA(2), 1)];
    [num2str(A(3), 1), ' \pm ', num2str(epsA(3), 1)];
    [num2str(VE*100, 3), '%']};

Trnsp = table(Quantity, Yuc_Sec, Flo_Sec)

%% Make the plots of the Figure 7

Ty = Tiy;
Py = Piy;

Tf = Tif;
Pf = Pif;

close all
figure; set(gcf, 'pos', [100 50.7778 1.6489e+03 696.8889], 'color', 'w');

subplot(3, 1, 1)
hold on;
plot(diy, smoothdata(Ty, 'gaussian', 7), 'linewidth', 1.3);

YYp = smoothdata(Ty(:), 'gaussian', 7) + smoothdata(Py(:), 'gaussian', 7);
YYm = smoothdata(Ty(:), 'gaussian', 7) - smoothdata(Py(:), 'gaussian', 7);
YY = [YYp(:); flipud(YYm(:)); YYp(1)];
XX = [diy(:); flipud(diy(:)); diy(1)];

hold on;
pt = patch(XX, YY, rgb('RoyalBlue'));pt.FaceAlpha = 0.4;
pt.EdgeColor = 'none';
rf = refline(0, 27.4); rf.Color = 'k';
datetick('x', 'mm-yyyy')
xlim([tid(1), tid(end)])
text(735068, 12.68379, '\textbf{(A)}', 'Interpreter', 'latex','FontSize', 15);
ylim([8.68, 41])
grid on; box on


subplot(3, 1, 2)
hold on;
plot(dif, smoothdata(Tf, 'gaussian', 7), 'r', 'linewidth', 1.3);

% create the uncertainty shade
YYp = smoothdata(Tf(:), 'gaussian', 7) + smoothdata(Pf(:), 'gaussian', 7);
YYm = smoothdata(Tf(:), 'gaussian', 7) - smoothdata(Pf(:), 'gaussian', 7);
YY = [YYp(:); flipud(YYm(:)); YYp(1)];
XX = [dif(:); flipud(dif(:)); dif(1)];

hold on;
pt = patch(XX, YY, rgb('Orange'));pt.FaceAlpha = 0.5;
pt.EdgeColor = 'none';
pt.EdgeColor = 'none';
rf = refline(0, 27.4); rf.Color = 'k';
datetick('x', 'mm-yyyy')
xlim([tid(1), tid(end)])
grid on; box on

text(735068, 12.68379, '\textbf{(B)}', 'Interpreter', 'latex','FontSize', 15);
ylim([8.68, 41])

subplot(3, 1, 3)
plot(diy, smoothdata(Py, 'gaussian', 7), 'color', rgb('RoyalBlue'), ...
    'linewidth', 1.3);
hold on;
plot(dif, smoothdata(Pf, 'gaussian', 7), 'color', rgb('red'), ...
    'linewidth', 1.3)
datetick('x', 'mm-yyyy')
xlim([tid(1), tid(end)])
hold on;
ylim([0, 5.5])
box on;
grid on;

ylabel('Error std (Sv)','Interpreter', 'latex','FontSize', 15)
text(735068, 1.31568, '\textbf{(C)}', 'Interpreter', 'latex','FontSize', 15);

