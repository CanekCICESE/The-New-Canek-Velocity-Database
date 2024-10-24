% funcion para preparar graficas de secciones vertiales de cualquier
% variable temperatura corriente promedios, desviacion stdr, modos, etc
% X and Z are the grid coordinates, da is the delta_area or bin size in
% squared meters

function [mask, xx, yy, xmap, zmap, celda, X, Z, da] = prep_sec(flg, vres, hres)
ir = pwd;
str = strfind(ir, '\');
if ~isempty(str)

    fu = 'prep_sec.m';
    a = which(fu);
    batroot = [a(1:strfind(a, fu)-1), '\Batis\'];

else
    batroot = [a(1:strfind(a, fu)-1), '/Batis/'];
end

switch flg
    case 'yuc'
        eval( [ 'load ', batroot, 'yucsec3.mat' ] );
        
        xlon=-86.8:hres:-84.9;
        profs=10:vres:2100; profs = -profs;
        [X,Z]=meshgrid(xlon,profs); 
        xx = xlon; yy = profs;
        
        xv = bp(:,1); yv = bp(:,3); xv = [xv ; xv(1)]; yv = [yv ; 0];
        
        celda = find(inpolygon(X(:),Z(:),xv,yv));
        mask=nan*ones(size(X)); dum=mask(:); dum(celda)=ones(size(celda));
        mask(:)=dum;
        
        celda=find(~isnan(mask));
        xmap=X(celda); zmap=Z(celda);
        
        dx = hres*3199.85/0.03;
        dz=vres; da=dx*dz;
        % area total de la seccion;
        
    case 'flo'
        eval( [ 'load ', batroot, 'flosec1.mat' ] );
        
        ylat=23.15:hres:24.65;
        profs=10:vres:1400; profs=-profs;
        nxb=length(ylat);
        
        % carga seccion del canal para graficas
        [X,Z]=meshgrid(ylat,profs);
        xv=bp(:,2); yv=bp(:,3); xv = [xv ; xv(1)]; yv = [yv ; 0];
       
        celda = find(inpolygon(X(:),Z(:),xv,yv));
        mask=nan(size(X)); dum=mask(:); dum(celda)=ones(size(celda));
        mask(:)=dum;
        xmap=X(celda); zmap=Z(celda);
        
        xx = ylat; yy = profs;
%         if gf
%             hold(ax, 'on');
%             plot(bp(:, 2),bp(:,3), 'color', rgb('Silver'),'linewidth',3);
%             axis([23 24.8 -1400 100]);
%         end
        iori=1; orisec=13.7825;  cosori=cosd(orisec); 
        %calcula parametros para transporte
        f1=50*60/27; f2=1000; % para convertir km a m
        facesc = f1*f2/cosori;
        dx=hres*f1*f2; if iori, dx=dx/cosori; end
        dz=vres; da=dx*dz;
           
end
% if ~gf
%     close all;
% end

end

