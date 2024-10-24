function [GL, ED, flgs] = GEL( x1, x2 )
% Obtain the effective degrees of freedom for two series s1 and s2
% considering the two criteria T_0 and T_e (e-folding) for the time decorrelation scales:
% ED = [T_0(s1), T_0(s2)
%      T_e(s1), T_e(s2)]
% in case one of the crossings does not exist, the point where it stabilizes is calculated, i.e., where dR/dx is approximately 0; therefore, few lags are used.
% The effective degrees of freedom are calculated using the equation from Panofsky and Brier (1958).
%
% GL are the GEL:
% GL = [GEL_0(s1), GEL_0(s2)
%      GEL_e(s1), GEL_e(s2)]
%
% flgs are the flags indicating whether there was (1) or was not (0) crossings of r through 0
% and exp(-1), respectively.
%
% Giovanni Durante, March 2021

% verifica el tamaño de las series
if length(x1) ~= length(x2)
   error('Las series deben tener el mismo número de datos') 
end

Nt = length(x1); % obtiene longitud de las series 
    
lgs = ceil(Nt/4); % un quinceabo de las series parece funcionar bien para evitar fuctuaciones espurio

[r1, ~] = autocorr(x1, lgs);
[r2, lags] = autocorr(x2, lgs);
R = [r1, r2];
flgs = [1, 1; 1, 1]; % banderas que indican cuando existe cruces por 0 y exp(-1)
dx = mean(diff(lags))/2;
for i = 1 : size(R, 2)
    
    r = R(:, i);
    T0 = polyxpoly( lags, r, [0, lags(end)], [0, 0] );

    if isempty(T0) % si no hay cruce por 0 toca calcular el punto donde se estabiliza derivando
        
        x = lags(2:end);
        y = r(2:end);
        
        xd = x(1:end-1)+dx;
        dy = abs(diff(y));
        [~, ind] = min(dy);
        T0 = xd(ind);
        T_0(i) = T0(1);
%         ['No existe cruce por cero, este es el punto donde se "estabiliza" T_0 = ', num2str(T0(1))]
        flgs(1, i) = 0;
    else % si hay cruce por 0, todo bien
        T_0(i) = T0(1);
%         ['Existe cruce por cero!. T_0 = ', num2str(T_0(i))]
        
    end
    % cruce por efolding
    Tef = polyxpoly( lags, r, [0, lags(end)], [exp(-1), exp(-1)] );
    if isempty(Tef) % si no hay
        x = lags(2:end);
        y = r(2:end);
        
        xd = x(1:end-1)+dx;
        dy = abs(diff(y));
        [~, ind] = min(dy);
        T0 = xd(ind);
        T_e(i) = T0(i);
%         ['No existe cruce por exp(-1), este es el punto donde se "estabiliza" T_0 = ', num2str(T_e(i))]
        flgs(2, i) = 0;
    else
        T_e(i) = Tef(1);
%         ['Existe cruce por exp(-1)!. T_e = ', num2str(T_e(i))]
    end
   ED(:, i) = [T_0(i); T_e(i)];
  
end

    GL(2, :) = (Nt*(dx*2))./(2*ED(2, :));
    GL(1, :) = (Nt*(dx*2))./(ED(1, :));
%     GL = (Nt*(dx*2))./(ED);

