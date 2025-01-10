%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
%
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% "MOLARVOLUME" solves the cubic EoS for Z (and V)
%==========================================================================
%
% INPUTS
% -- T = temperature in K [REAL 1x1]
% -- P = pressure in bar [REAL 1x1]
% -- b = co-volume, L/mol [REAL NCx1]
% -- abRT = a/(b*R*T) of the mixture [REAL 1x1]
% -- c = volumetric translation, L/mol [REAL NCx1]
% -- phase = 1 (liquid) or 2 (vapor) [INTEGER 1x1]
%
% OUTPUTS
% -- Z_res = Compressibility factor of the mixture [REAL 1x1]
% -- V_res = molar volume of the mixture in L/mol [REAL 1x1]
%
%==========================================================================

function [Z_res,V_res] = MOLAR_VOLUME(T,P,b,abRT,phase)

%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2;

% atractive parameter
a = (b*R*T)*abRT;

% dimensioless EoS parameters 
A = a.*P./(R.*T)^2;
B = b.*P./(R*T);

% Compressibility factor
Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
ZR = [];
for i = 1:3
   if isreal(Z(i)) && Z(i) > 0
   	    ZR = [ZR real(Z(i))];
   end
end

% Picking the right root...
if phase == 1 % (liquid)
    Z_res = (min(ZR));

elseif phase == 2 % (vapor)
    Z_res = (max(ZR));

else % most stable =  lower fugacity
    lnphi = (ZR - 1) - log(ZR - B) - A/(2 * B * sqrt(2)) * log((ZR + (1 + sqrt(2)) * B)./(ZR + (1 - sqrt(2)) * B));
    [~,idx] = min(exp(lnphi));
    Z_res = ZR(idx);
end

% Molar volume for the chosen phase
V_res = (R*T/P)*Z_res;

%==========================================================================
end