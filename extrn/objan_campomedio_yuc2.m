 function [maps, Wlarge, Err]=objan_campomedio_yuc2(data,xdat,zdat,params,xmap,zmap)
%%%%function [maps,Wlarge,datamean]=objan_campomedio_yuc(data,xdat,zdat,params,xmap,zmap);
%   function maps=objmaps(data,xdat,zdat,params,xmap,zmap);
%                                               J. Ochoa , Jul, 99
%      params=[largex largez smallx smallz largevar smallvar];
% Pp
%            objective mapping function
% Pp
% on exit there are as many maps as columns of "data",
% "xdat" & "zdat" are the position (columns) coordinates of "data", 
% the "maps" will be estimated values on the gridd defined by 
% "xmap" & "zmap"  and will depend on the SIX parameters given
% in array "params"=[largex largez smallx smallz largevar smallvar];
% w/ large.. parameters => an estimation of the data and
% corresponding residuals is done on the positions of "data"
% (i.e. on (xdat,zdat)), then
% w/ small.. parameters a refinement is done which includes small 
% scale variations.
%first column of maps is the map, second column is the error
%the error is missing a term from the estimation of the mean.
%JS modificado de codigo de Pepe
%Se aniade matrz de salida covarerr que es la matriz de covariancia del
%error de estimacion se aniadiero calculos al final
%JS aug 2001
% 
% se modifica para ajustar solo el campo medio
% JC 05/2016


f1=50*60/27; centlat=21.63;
f2=cosd(centlat)*1000;
% factor de escala por la orientacion de la seccion, por consistencia
orisec=10.9096;  cosori=cosd(orisec); 
facesc = f1*f2/cosori;

% compute correlation matrices based on positions and parameters
m=length(xdat);
xx=xdat';xx=xx(ones(m,1),:)-xdat(:,ones(1,m));xx=(xx*facesc)/params(1);xx=xx.*xx;
zz=zdat';zz=zz(ones(m,1),:)-zdat(:,ones(1,m));zz=zz/params(2);zz=zz.*zz;
ALAsru=exp(-(xx+zz)); ALAinv=ALAsru + diag(params(3)*ones(1,m)); ALAinv=inv(ALAinv); 

% ajustar un plano a los datos para eliminarlo antes de la IO
%x=facesc*xdat(:)/params(1); z=zdat(:)/params(2);
%A=[ones(m,1),x,z];
%A=[ones(m,1),x,z,x.^2,z.^2,x.*z];
%c=A\data, dataplano=A*c; nx=length(xmap); 
%xm=facesc*xmap(:)/params(1); zm=zmap(:)/params(2);
%plano=[ones(nx,1),xm,zm]*c;
%plano=[ones(nx,1),xm,zm,xm.^2,zm.^2,xm.*zm]*c;
%data=data-dataplano;
% compute  means & substract from data
%datamean=sum(ALAinv*data)/sum(sum(ALAinv)),
%%%%datamean=mean(data); data=data-datamean;

% compute residuals & weights
Wlarge=ALAinv*data; %size(ALAinv); %size(data), size(Wlarge), size(ALAsru),

xmap=xmap(:); zmap=zmap(:); nx=length(xmap);
xp=xmap'; zp=zmap'; xx=xdat; zz=zdat;
px=xx(:,ones(1,nx))-xp(ones(m,1),:); 
pz=zz(:,ones(1,nx))-zp(ones(m,1),:); 
px1=(facesc*px/params(1)).^2; pz1=(pz/params(2)).^2; %max(max(abs(px1))),max(max(abs(pz1)))
expdl=exp(- (px1+pz1)); 
%%%%maps(:,1)= expdl'*Wlarge + datamean; 
maps(:,1)= expdl'*Wlarge;

%Computing normalized squared percentage error relative to the "true"
%field x

Err = nan(size(xmap));
numErr = expdl'*ALAinv*expdl;
% Err(:,1)= expdl'*Ainv*expdl;
for i=1:length(numErr)
      Err(i)=eye(length(numErr(i,i)))-numErr(i,i);
end






%maps(:,1)= expdl'*Wlarge + plano; size(expdl)

