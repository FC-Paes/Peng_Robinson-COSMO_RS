%==========================================================================
% "MIX_RULES" calculates the mixture parameters of the cubic EoS
%==========================================================================
%
% INPUTS
% -- z = vector containing the mole fractions of the Ncomp components [REAL NCx1]
% -- b = vector of pure compounds co-volumes, L/mol [REAL NCx1]
% -- c = vector of pure compounds volumetric translation, L/mol, L/mol [REAL NCx1]
% -- abRT = vector of pure compounds "a/(bRT)" values [REAL NCx1]
% -- lnGam = vector of ln(activity coefficient of i) [REAL NCx1]
%
% OUTPUTS
% -- b_mix = co-volume of the mixture in L/mol [REAL 1x1]
% -- abRT_mix = a/(RTb) of the mixture, dimensionless [REAL 1x1]
% -- der_abRT_mix =  vector containing the derivatives of abRT_mix wrt. z(i) [REAL NCx1]
% -- sum_zb vector = vectot containing the linear forms associated to bmix
%                    ie, sum over j of z(j) * bij(i,j) in L/mol [REAL NCx1]
%
%==========================================================================

function [b_mix,c_mix,abRT_mix,sum_zb,der_abRT_mix] = MIX_RULE_EOS(z,b,c,abRT,lnGam,par)

% number of compounds
NC = length(z);

% mix rule constant
q1 = par.univ_mixrule.q1;

%--------------------------------------------------------------------------
% Calculation of the mixture co-volume
%--------------------------------------------------------------------------

bij = zeros(NC,NC);
% bij matrix
s = par.univ_mixrule.s;
for i = 1:NC
    for j = i:NC
        bij(i,j) = ((b(i)^(1/s) + b(j)^(1/s))/2)^s;
        if (i~=j)
            bij(j,i) = bij(i,j);
        end
    end
end
% b_mix
b_mix = dot(z,bij*z);

% linear form associated with the quadratic function bmix
sum_zb = (z'*bij)';

%--------------------------------------------------------------------------
% Calculation of the mixture volume translation
%--------------------------------------------------------------------------
c_mix = dot(z,c);

%--------------------------------------------------------------------------
% Calculation of the mixture attractive parameter
%--------------------------------------------------------------------------

% excess Gibbs energy
Ge_RT = dot(z,lnGam);

% sum of z*abRT
sum_abRT = dot(z',abRT);

% abRT_mix
abRT_mix = Ge_RT/q1 + sum_abRT;

% derivative of abRT_mix wrt. z(i)
der_abRT_mix = abRT + lnGam./q1;

%==========================================================================
end