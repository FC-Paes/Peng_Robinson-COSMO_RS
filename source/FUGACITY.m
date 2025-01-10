%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
%
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% "FUGACITY" calculates the ln of the fugacity coeffcient of each compound
% in the mixture
%==========================================================================
%
% INPUTS
% -- T = temperature in K [REAL 1x1]
% -- P = pressure in Pa [REAL 1x1]
% -- V = molar volume in L/mol [REAL 1x1] (already corrected with the volumetric  translation)
% -- b_mix = co-volume of the mixture in L/mol [REAL 1x1]
% -- abRT_mix = a/(RTb) of the mixture, dimensionless [REAL 1x1]
% -- der_abRT_mix =  vector containing the derivatives of abRT_mix wrt. z(i) [REAL NCx1]
% -- sum_zb vector = vectot containing the linear forms associated to bmix
%                    ie, sum over j of z(j) * bij(i,j) in L/mol [REAL NCx1]
% -- c = pure compound volumetric translation [REAL NCx1]
%
% OUTPUTS
% -- FUG =ln(fugacity coeff) of each compound in mixture [REAL NCx1]
%
%==========================================================================

function FUG = FUGACITY(T,P,V,b_mix,sum_zb,der_abRT_mix,c)
%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2;
% ln(fugacity coefficient) = FUG
F1 = log(R.*T./(P.*(V - b_mix)));
F2 = (2.*sum_zb./b_mix - 1).*(P.*V/(R.*T) - 1);
sqrt2 = sqrt(2);
F3 = - der_abRT_mix./(2.*sqrt2).*log((V + (1 + sqrt2).*b_mix)./...
    (V + (1 - sqrt2).*b_mix));
FUG = F1 + F2 + F3;
FUG = FUG + P.*c./(R.*T);

%==========================================================================
end