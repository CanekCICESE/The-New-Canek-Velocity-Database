% get the bathymetric profile of the Canek sections to plot it  
% inputs: 
% flg: 'yuc' or 'flo' determines which section we want to plot
% ax: is the axis, normally gca
% color: rgb triplete or matlab colorcode 
% lw: linewidth of the profile

function [xto, yto] = batsec(flg, ax, color, lw)
    
    A = pwd;
    fu = 'prep_sec.m';
    a = which(fu);
    
    try
        batroot = [a(1:strfind(a, fu)-1), '\Batis\'];
        cd(batroot);
    catch
        batroot = [a(1:strfind(a, fu)-1), '/Batis/'];
    end
    cd(A)
    switch flg
        case 'yuc'
            load([batroot, 'yucsec3.mat']);
            xto = bp(:,1); yto = bp(:,3);
            hold(ax, 'on');
            plot(xto, yto,'color', color,'linewidth', lw);
            axis([-86.8 -84.9 -2100 100]);

        case 'flo'
            load([batroot, 'flosec1.mat']);
            xto = bp(:,2); yto = bp(:,3);
            hold(ax, 'on');
            plot(xto, yto,'color', color,'linewidth', lw);
            axis([23 24.8 -1400 100]);           
    end
end

