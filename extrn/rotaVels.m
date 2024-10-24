function [cn, cp] = rotaVels(flag, u, v)

switch flag

    case 'yuc'
        orisec = 10.9096; % 9.3270; 10.9096;
        un = cosd(orisec + 90); vn = sind(orisec + 90); % normal vector to the section
        up = cosd(orisec); vp = sind(orisec); % parallel vector to the section

    case 'flo'
        orisec = 15.0511;
        up = cosd(orisec + 90); vp = sind(orisec + 90); % normal vector to the section
        un = -cosd(orisec + 180); vn = -sind(orisec + 180); % parallel vector to the section
end

cn = nan([size(u)]);
cp = nan([size(u)]);

for kk = 1 : size(u, 2)
    for jj = 1 : size(u, 1)
        cn(jj, kk) = [u(jj, kk), v(jj, kk)] * [un, vn]'; % normal component
        cp(jj, kk) = [u(jj, kk), v(jj, kk)] * [up, vp]'; % parallel component
    end
end
end
