% Example of the observational noise std estimation
% this code uses the  function mspec from  the Jlab library:
% You can download the Jlab from: https://jmlilly.net/code.html or https://github.com/jonathanlilly/jLab

close all; clc; clear all;

% path to the file
yourPath = 'C:/Users/nucle/Tesis/Papers/Paper_I/Figures/Programas/CreaRepos/';
origin = [yourPath, '/Data/BackgroundNoiseSTD/'];

seccion = 'Yucatan';
% seccion = 'Florida';

flag = lower(seccion);
flag = flag(1:3);

% possible moorings
switch seccion
    case 'Yucatan'
        ancla = {'YUC1', 'YUC2', 'YUC3', 'YUC4', 'YUC5', ...
            'YUCI5', 'YUC6', 'YUCI6', 'YUC7', 'YUCI7', 'YUC8', ...
            'YUCI8', 'YUC9', 'YUCI9', 'YUC10'};
        ani = 1:length(ancla);
    case 'Florida'
        ancla = {'EFL1', 'EFL2', 'EFL3', 'EFL4', 'EFL5', 'EFLI5', 'EFL6', 'EFL7'};
end

% canek campaigns
CNK = [26,  29,  34, 37,  39,	42,  48];


an = 9; % mooring, for this example is YUC7
kk = 3; % canek campaign, for this exam,ple is CNK34
cnk = num2str(CNK(kk));

anc = ['CNK', cnk, '_', char(ancla(an))];
load([origin, anc])

cou = 1;
% anchor is an array of structures. hence, the i-est element of anchor: 'anchor(i)'
% is a structure containing all the data and information from a given
% instrument in the mooring 'an' during the campaing 'kk'

% is anchor exists and is not empty, proceed:
if exist('anchor') && ~isempty(fieldnames(anchor))

    for lu = 1 : 1 length(anchor)
        %
        t = anchor(1).time;
        name = anchor(lu).name;
        tip = find(int16(name) == 45);
        rtipo = char(name( tip(2)+1:tip(3)-1 ));
        z = -abs(anchor(lu).z);
        p = abs(anchor(lu).p);
        nm = size(z);

        T = anchor(lu).T;
        X = anchor(lu).x;
        Y = anchor(lu).y;
        V = anchor(lu).v;
        U = anchor(lu).u;
        W = anchor(lu).w;

        if nm(2) > nm(1)
            V = nanmean(V, 1);
            U = nanmean(U, 1);
            W = nanmean(W, 1);

        else
            V = nanmean(V, 2);
            U = nanmean(U, 2);
            W = nanmean(W, 2);
        end

        X = nanmean(X(:)); Y = nanmean(Y(:));

        inba = isnan(U(:) + V(:));
        U(inba) = []; V(inba) = []; t(inba) = [];
        U = detrend(U(:), 1); V = detrend(V(:), 1);
        %
        if ~isempty(V)
            clear i;

            dt = 60*60;
            cv = V;
            psi=sleptap(length(cv),5);

            [f,spp]=mspec(1,cv,psi,'cyclic'); % we use Jlab fncs to compute the power spectral density estimation
            df = mean(diff(f));
            %Compute the power spectral density estimation:
            [fd,sppd]=mspec(1/24,cv,psi,'cyclic');

            xo = flipud(f(:));
            xf = flipud(fd(:));
            ps = flipud(spp(:));

            li = find( xf >= 2.5);

            xo = xo(li);
            xf = xf(li);
            ps = ps(li);

            psf = detrend(ps(li),0);
            psf = psf./max(psf);
            xff = (xf);

            % fit a general gaussian curve
            fi = fittype('a*(exp((x.^2).*b)) + c');
            [fite,gof,fitinfo] = fit(xff, psf, fi, 'StartPoint',[1, 0, 0]);
            fua = feval(fite, xf);
            L = coeffvalues(fite);
            d = L(3); % vertical shift part

            thresh = 0.02; % threshold to look for the flat part
            ini = find(fua-d <= thresh );


            % Using the trapezoidal integration to find the variance
            % estimated noise variance integrating the flat part
            std_noise(cou) = sqrt(abs(trapz(ps(ini), xo(ini))));
            inis = (fua-d <= thresh );
            
            % estimated signal variance total - noise variance
            var_to_noise(cou) = ( (abs(trapz(ps(ini), xo(ini)))))/ (var(cv) - (abs(trapz(ps(ini), xo(ini)))));

            % for graphic porpuses =========================================================================
            psf = ps(li);
            xff = xo(li);
            sill = max(psf);
            
            % fit a general gaussian curve
            fi = fittype('a*(exp((x.^2).*b)) + c');
            [fite,gof,fitinfo] = fit(xff, psf./sill, fi, 'StartPoint',[max(psf), 0, 0]);
            fua = feval(fite, xff);


            figure; set(gcf, 'Position', ...
                [82.7778 27.6667 682.6667 832])

            subplot(2, 1, 2)
            plot(xff, psf*1e4, 'linewidth', 1.1); hold on;
            plot(xff(ini), psf(ini)*1e4, 'k', 'linewidth', 1.1); hold on;
            grid on;
            cax= ylim;

            plot(xff, fua*sill*1e4, 'LineWidth', 2, 'color', rgb('OrangeRed'));
            xlabel('Frequency (hrs$^{-1}$)', 'interpreter',  'latex');
            ylabel('Power spactral density (m s$^{-1}$ hrs)', 'interpreter',  'latex');
            set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 14)
            text(0.458712, 1.36720, '\textbf{(II)}', 'interpreter',  'latex', 'fontsize', 18);

            text(0.26750, 1.0262, ...
                ['Estimated error STD: ', num2str(std_noise*100, 3), ' cm s$^{-1}$'], 'interpreter',  'latex', 'fontsize', 15);

            subplot(2, 1, 1)
            set(gca, 'pos', [0.13 0.519673711203466 0.775 0.341162790697675])
            loglog(f,spp, 'linewidth', 1.1); hold on;
            loglog(xo(ini), ps(ini), 'k', 'linewidth', 1.1);
            plot(xff, fua*sill, 'LineWidth', 3, 'color', rgb('OrangeRed'));
            grid on;
            xlabel('Frequency (hrs$^{-1}$)', 'interpreter',  'latex');
            ylabel('Power spactral density (m s$^{-1}$ hrs)', 'interpreter',  'latex');
            set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 14)
            xlim([0.000132214079303354 0.661070396516769])
            text(0.249152377523668, 0.437799946244678, '\textbf{(I)}', 'interpreter',  'latex', 'fontsize', 18);

            set(gca, 'YTickLabelRotation', 50)

            annotation('line', 'Position',[0.734126984126984 0.625475285171103, ...
                -0.600793650793653 -0.17490494296578])

            annotation('line', 'Position',[0.882539682539682 0.583650190114068,...
                0.0180555555555557 -0.13296102661597])

            title({['NS10933 LongRanger ADCP in YUC7 CNK34 at $z$ $\approx$ 170 m']}, 'interpreter',  'latex');

            ll = legend({'Power spectral density', 'Flat segment (white noise)', 'Gaussian fit'}, 'interpreter',  'latex');
            ll.Position = [0.141395494595641 0.629596004014686 0.349480417813118 0.0700237678484316];


        end
    end

end
%

%
