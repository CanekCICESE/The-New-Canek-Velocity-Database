function [A,Fase,epsA,epsFase,dta,VE,as,epsas]=ajustew(t,d,Narm,T,W) 
% Harmonic adjustmen
%------------------------------------------------------------------------
%  The output is vectors    A(1) mean value
%                         Fase(1) phase of the mean value
%                            A(2) amplitude of the first harmonic
%                         Fase(2) phase of the first harmonic, etc. etc......
%                                                          Reieb El Cotur
%------------------------------------------------------------------------
% AJUSTEW Harmonic fit: A cos(wt - Phase) 
%     [A F epsA epsF dta VE a epsa]=AJUSTEW(t,d,Narm,T,[W])
%     INPUT:
%                     t - times
%                     d - data 
%                  Narm - number of harmonics
%                     T - period
%                    [W - symmetric and positive definite weight matrix
%                         (Default: W = I)]
%                  nargin= < 5 no weights
%     OUTPUT:
%                A,espA - amplitude and uncertainty
%          Fase,epsFase - phase and error
%                   dta - residue/(standard deviation data)
%                    VE - explained variance (= 1 - <dta^2>)
%                a,epsa - coefficients and uncertainties d=a*F+r
%     It also works if d is a matrix array.

%------------------------------------------------------------------------
% Arrange t and d to fit d=as*F+r
nargin
[Nt nt]=size(t); 
if Nt~=1 t=t'; nt=Nt; end
[N n]=size(d); 
if n~=nt; d=d'; [N n]=size(d); end
m=2*Narm+1;                           % number of functions  
wt=2*pi/T*[1:Narm]'*t;                % phases
F=[ones(size(t));cos(wt);sin(wt)];    % functions
%
% fitting and error estimation
if nargin<5 % without weights
  as=d/F;                             % estimated coefficients  
  r=d-as*F;                           % residues
  sg2=...                             % squared error per point 
    diag(r*r')...                     % sum of squared residues
    /(n-m);                           % degrees of freedom
  Ca=diag(inv(F*F'))*sg2';            % co-errors coefficients
  epsas=sqrt(Ca);                     % estimated errors 
  dm=mean(d')'*ones(1,n);             % mean data
  sga=sqrt(mean((d-dm)'.^2));         % standard deviation 
  dta=diag(1./sga)*r;                         
  VE=(1-mean(dta'.^2))';              % explained variance
else        % with weights
  G=F*sqrt(W);
  as=d*sqrt(W)/G;                       
  r=d-as*F;                           
  sg2=...                              
    diag(r*W*r')...                   % weighted sum of squared residues     
    /(n-m);                           
  Ca=diag(inv(G*G'))*sg2';            
  epsas=sqrt(Ca);                      
  dm=meanw(d,W)'*ones(1,n);           % weighted mean data
  sga=sqrt(meanw((d-dm).^2,W));       % weighted standard deviation
  dta=diag(1./sga)*r;                         
  VE=(1-meanw(dta.^2,W))';            
end

% amplitudes and phases
as=as'; 
a=as(1:Narm+1,:);              
b=[zeros(1,N); as(Narm+2:2*Narm+1,:)];
epsa=epsas(1:Narm+1,:);              
epsb=[zeros(1,N); epsas(Narm+2:2*Narm+1,:)];
[A Fase epsA epsFase]=afase(a,b,epsa,epsb);
end







function [A,F,eA,eF]=afase(a,b,ea,eb)
% AFASE Amplitude, phase, and uncertainties.
%   [A Fase epsA epsFase]=AFase(a,b,epsa,epsb)
%   Input: a, b, uncertainties    [a cos wt + b sin wt]
%   Output: A, Phase, uncertainties  [A cos(wt - Phase)]
%   Works with matrices.
A=sqrt(a.^2+b.^2);
F=atan2(b,a);
eA=sqrt((a.*ea).^2+(b.*eb).^2)./A;
eF=sqrt((b.*ea).^2+(a.*eb).^2)./A.^2;
end
