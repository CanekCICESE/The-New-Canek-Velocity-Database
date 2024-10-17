% Used to create the correlation functions show in
% the subsection 2c, Figure 3.
%
% G. Durante, 2024
clear; close all; clc;

% path to the main folder of the repository
raiz = 'C:\Users\nucle\Tesis\Papers\Paper_I\Figures\Programas\CreaRepos\';

CNKi = [29  34	37  39	42  48]; % canek campaigns (2012-2020)
filtra = 0;

ch = [];
hs = [];
id = [];
flag = 'yuc';



%% Compute horizontal correlogram

for ll = 1 : length(CNKi)

    load(...
        [raiz, 'Data\CampaignAverages\',...
        'Clim_yuc_vel_', num2str( CNKi(ll) ), '_raw'])

    yp = getCanekSCoor(flag, xp); % obtain the latitudes of the observations
    xp = xp(:); yp = yp(:); zp = zp(:); vn = vn(:);

    [xp, aii] = sort( xp, 'ascend' );
    yp = yp(aii); vn = vn(aii); zp = zp(aii);

    %=========================================================
    % some times, filtering the profiles average a bit is appropriate
    [xu, aii] = unique(xp);
    yu = yp(aii);

    if filtra
        for k = 1 : length(xu)
            indp = find(xp == xu(k));
            vn(indp) = smoothdata(vn(indp), 'gaussian', 10);
        end
    end

    dx1 = 0.15;
    dx = [0:dx1:1.8];
    delz = 0;

    vn = vn - mean(vn);
    xco = xp; zco = zp; vco = vn;

    close all;
    figure; set(gcf, 'pos', [50 50 1536 748.8]);
    axes('pos', [0.0609375 0.583837209302326 0.234798180006381 0.341162790697675])
    scatter(xco, zco, [], vco, 'filled')
    hold on

    axis([-86.3296767116409 -85.3912904035724, -1030.9375 55.3125])
    box on;
    grid on;
    cb = colorbar();
    cb.Position = [0.299402845131631 0.583333333333333 0.0138888888888889 0.341880341880342];

    plot(-86.229, -400, 'sk', 'markerfacecolor', 'k');

    plot(-85.56998, -400, 'xk', 'linewidth', 1.2);
    plot(-86.229+(dx(5+1)), -400, '>k', 'markerfacecolor', 'k');
    quiver(-86.229, -400, dx(5+1)+0.08, 0, 'k');


    hold on;
    dve = [-86.229, -86.229+0.16];
    diff(dve)
    quiver(-86.229, -254.48, 0.16, 0, 'k');

    ax0(1) = plot(-86.229, -254.48, 'sk', 'markerfacecolor', 'k');
    ax0(2) = plot(-86.14104, -254.48, 'xk', 'linewidth', 1.2);
    ax0(3) = plot(-86.229+dx(5), -400, '*k', 'linewidth', 1.2);
    ax0(4) = plot(-86.229+0.15, -254.48, '>k', 'markerfacecolor', 'k');

    texto = text(-86.179, -175.518, '$r^1_x$', 'interpreter', 'latex', 'fontsize', 15);
    texto = text(-85.918, -338.448, '$r^5_x$', 'interpreter', 'latex', 'fontsize', 15);

    leg = legend(ax0, {'Tail values', 'Head values', 'Lower bin bound', 'Upper bin bound'}, 'interpreter','latex');
    leg.Position = [0.18266621009313 0.815482547565063 0.106743885839957 0.103498933916418];
    title('\textbf{(A)}', 'interpreter','latex');
    ylabel(cb, '(m/s)', 'interpreter','latex')


    rhox = [];
    % ch = [];
    % hs = [];

    for j = 1 : length(dx)-1
        px = [];
        py = [];
        pz = [];
        pl = [];
        for i = 1: length(vco)
            lt = xco; zt = zco; tet = vco;
            %         pivot.XData = lt(i); pivot.YData = zt(i);
            lt(i) = []; zt(i) = []; tet(i) = []; % quitar el punto que tomamos como pivote
            xclu = (lt - xco(i));
            %         xclu = xclu(xclu < 0);
            iex = find( xclu < dx(j+1) & xclu >= dx(j)+0.05  );
            lt = lt(iex);
            zt = zt(iex);
            tet = tet(iex);

            zclu = abs(-zt + zco(i));
            iex = find( zclu <= delz );
            lt = lt(iex);
            zt = zt(iex);
            tet = tet(iex);

            mn = length(lt);
            if mn > 0
                if  mn > 1
                    ptx = repmat(vco(i), [1, mn]);
                    ptl = repmat(xco(i), [1, mn]);
                    ptz = repmat(zco(i), [1, mn]);
                else
                    ptx = vco(i);
                    ptz = zco(i);
                    ptl = xco(i);
                end
                px = [px; ptx(:)];
                py = [py; tet(:)];
                pz = [pz; ptz(:)];
                pl = [pl; ptl(:)];
            end
            %         [i, length(xco), j, length(dx)-1] % print the progress


        end
        %========================================================================================
        % All this block is just to plot an example of the correlation at
        % different distances and the scatter plot between heads and tails

        if j == 1
            subplot(2, 3, 2)
            plot(px, py, 'ok', 'markerfacecolor', rgb('RoyalBlue'));
            axis([-0.4, 1, -0.4, 1])
            box on;
            grid on;
            co = num2str(corr(px, py), 1);
            title('\textbf{(B)}', 'interpreter','latex')
            xlabel('Tails (m/s)', 'interpreter','latex')
            ylabel('Heads (m/s)', 'interpreter','latex')

            texto = text(-0.36063,  0.915781, ...
                ['$r^1_x = (0, 15]$ km. Corr = ', co], 'interpreter', 'latex', 'fontsize', 15);
        end
        if j == 6
            axes('pos', [0.635604619565214 0.583837209302326 0.213405797101447 0.341162790697675])
            plot(px, py, 'ok', 'markerfacecolor', rgb('RoyalBlue'));
            axis([-0.4, 1, -0.4, 1])
            box on;
            grid on;
            clc
            co = num2str(corr(px, py), 1);
            title('\textbf{(C)}', 'interpreter','latex')
            xlabel('Tails (m/s)', 'interpreter','latex')
            texto = text(-0.36063,  0.915781, ...
                ['$r^6_x = (75, 80]$ km. Corr = ', co], 'interpreter', 'latex', 'fontsize', 15);
            set(gca, 'YtickLabels', {})
        end
        %========================================================================================
        xul = unique(pl);
        aaa = 1;
        for ju = 1 : length(xul)
            indl = find( pl == xul(ju) );
            if length(indl) > 8
                cori(aaa) = nanmean((py(indl) - px(indl)).^2); % mean squared differences only if there are more than 8 cases
                aaa = aaa + 1;
            end
        end
        rhox(j) = nanmean(cori); % structure function
    end
    % estimate the bin size in meters
    dh = spheric_dist(yu(1), yu(1), xu(1), xu(1)+1);
    hdist = dx*dh;
    %=======================================================
    % add zeros at dx = 0 for a better fit
    rhox = smooth([0; rhox(:)]); % using smooth to reduce spikes a bit in the structure function
    % compute the sill and apply the relationship between the structure function or variogram
    %  and the correlation. (Bohling 2005; Webster and Oliver 2007; Kleijnen 2017)
    sill = max(rhox); % estimate the sill
    Ch = (sill - rhox)/(sill); % Compute the correlogram from the variogram
    % concatenate qith the rest of Canek campaigns
    ch = [ch(:); Ch(:)]; % correlogram
    hs = [hs(:); hdist(:)]; % distances
    %=======================================================
end


set (findobj (gcf, '-property', 'interpreter'), 'interpreter','latex')
set (findobj (gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter','latex')

Fs = findall(gcf,'-property','FontSize');
for k = 1 : length(Fs)
    if strcmp(get(Fs(k), 'type'), 'axes')
        if Fs(k).FontSize < 12.5
            Fs(k).FontSize = 12.5;
        end
    end
end


%%

figure('pos', [50 50 1.1844e+03 488]);
subplot(1, 2, 1)
hold on;
plot(hs(:)./1000, ch(:), 'sk', 'markerfacecolor',...
    [0.5273, 0.8047, 0.9792], 'markersize', 10);
hold on;

% fit gaussian curve:
xi = [0:100:max(hdist)+40000]*1e-3;
f = fittype('(exp((x.^2).*b))');
[fite,gof,fitinfo] = fit(hs./1000,ch,f, 'StartPoint',0);

L = coeffvalues(fite);
L = sqrt( 1/abs(L) );
cova = feval(fite, xi);
hold on;
plot(xi, cova, 'Color', [0.2734, 0.5078, 0.7031], 'linewidth', 2);
plot(xi, exp( -( ( xi./70 ) .^2) ), 'k', 'linewidth', 2);

xlabel('Horizontal distance');
ylabel('correlation distance');
title('Horizontal correlogram');
ll = legend({'Empirical', 'LS fit', 'Volume-conserv.'});
% save corre_yuc hs ch

%% Compute vertical correlogram


cv = [];
hv = [];

for ll = 1:6

   load(...
        ['C:\Users\nucle\Tesis\Presentaciones\OSM_2024\Presentacion\Figuras\DataNMethods\',...
        'Proms_x_canek\Clim_yuc_vel_', num2str( CNKi(ll) ), '_raw'])

    
    yp = getCanekSCoor(flag, xp); % obtain the latitudes of the observations
    xp = xp(:); yp = yp(:); zp = zp(:); vn = vn(:);
    vn = vn - mean(vn);
    xco = xp; zco = zp; vco = vn;

    dzi = 40;
    dz = [0:dzi:850];
    delx = 0;
    pinta = 0;

    rhoz = [];
    yot = [];
    nnz = [];
    vardifez = [];
    Cv  = [];

    for j = 1 : length(dz)-1

        px = [ ];
        py = [ ];
        pz = [ ];
        pl = [ ];
        for i = 1: length(vco)

            lt = xco; zt = zco; tet = vco;
            lt(i) = []; zt(i) = []; tet(i) = []; % quitar el punto que tomamos como pivote

            xclu = (zt - zco(i));
            iex = find( xclu < -dz(j) & xclu >= -dz(j+1) );

            lt = lt(iex);
            zt = zt(iex);
            tet = tet(iex);

            xclu = abs(lt - xco(i));
            iex = find( xclu <= delx );

            lt = lt(iex);
            zt = zt(iex);
            tet = tet(iex);

            if pinta
                pt.XData = xco(i);
                pt.YData = zco(i);
                pt2.XData = lt;
                pt2.YData = zt;
                drawnow
            end

            mn = length(lt);
            if mn > 0
                if  mn > 1
                    ptx = repmat(vco(i), [1, mn]);
                    ptl = repmat(xco(i), [1, mn]);
                    ptz = repmat(zco(i), [1, mn]);
                else
                    ptx = vco(i);
                    ptz = zco(i);
                    ptl = xco(i);
                end
                px = [px; ptx(:)];
                py = [py; tet(:)];
                pz = [pz; ptz(:)];
                pl = [pl; ptl(:)];
            end
            [i, length(xco), j, length(dz)-1]
        end

        xul = unique(pz);
        aaa = 1;
        for ju = 1 : length(xul)
            indl = find( pz == xul(ju) );
            if length(indl) > 6
                % cori(aaa) =  var( px(indl) - py(indl) );
                cori(aaa) = nanmean((py(indl) - px(indl)).^2);
                aaa = aaa + 1;
            end
        end
        rhoz(j) =  nanmean(cori);
    end

    vdist = dz(:)+20;
    rhoz = [0; rhoz(:)];

    sill = max(smooth(rhoz));
    Cv = (sill - smooth(rhoz))/(sill);
    vdist(Cv == 0) = [];
    Cv(Cv == 0) = [];

    cv = [cv(:); Cv(:)];
    hv = [hv(:); vdist(:)];
end

%%
subplot(1, 2, 2)

hold on;
plot(hv(:), cv(:), 'sk', 'markerfacecolor',...
    [0.5273 0.8047 0.9792], 'markersize', 10);
hold on;

%
xi = [0:1200];

f = fittype('(exp((x.^2).*b))');
[fite,gof,fitinfo] = fit(hv,cv,f, 'StartPoint',0);

L = coeffvalues(fite);
L = sqrt( 1/abs(L) );
cova = feval(fite, xi);
hold on;
plot(xi, cova, 'Color', [0.2734 0.5078 0.7031], 'linewidth', 2);
plot(xi, exp( -( ( xi./500 ) .^2) ), 'k', 'linewidth', 2);
%%
% save corre_yuc hs ch hv cv

