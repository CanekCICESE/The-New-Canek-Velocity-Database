% obtaine the latitude or longitude of the Canek moorings given the initial
% coordinate. 
% G. Durante, 2024

function ocoor = getCanekSCoor(flag, icoor)
switch flag
    case 'yuc'
        B = [0.1642   35.7545];
        ocoor = polyval(B, icoor);
    case 'flo'
        B = [-0.2689  -74.8350];
        ocoor = polyval(B, icoor);
end