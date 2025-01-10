%==========================================================================
% "PARAMETERS_EoS" calculates the parameters a(T) and b of the cubic EoS
% for pure compounds
%==========================================================================
%
% INPUTS
% -- T = temperature in K [REAL 1x1]
% -- P = pressure in bar [REAL 1x1]
% -- Tc = critical temperature in K [REAL NCx1]
% -- Pc = critical pressure in K [REAL NCx1]
% -- w = accentric factor [REAL NCx1]
% -- Twu 91 parameters (L,M, and N) [REAL NCx1]
% -- c = pure compound volumetric translation (if available) [REAL NCx1]
%                                             (if it is not available, c = 8888)
%
% OUTPUTS
% -- a = attractive parameter in ... [REAL NCx1]
% -- b = co-volume in m^3 [REAL NCx1]
% -- abRT = a/(bRT) [REAL NCx1]
% -- c = pure compound volumetric translation [REAL NCx1]
% (equals to input if available or calculated by correlation if not available)
%
%==========================================================================

function [a,b,c,abRT,alpha] = PARAMETERS_EoS(T,Tc,Pc,w,L,M,N,c,alpha_function)

%--------------------------------------------------------------------------
% Number of compounds
NC = size(Tc,1);

%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2; 

% Reduced temperature
Tr = T./Tc;

%--------------------------------------------------------------------------
% Parameters of the EOS for a pure component
%--------------------------------------------------------------------------

% Alpha function
switch alpha_function
    case 'TWU91'
        alpha = ones(NC,1);
        for i=1:NC
            if L(i,1) == 8888 || M(i,1) == 8888 || N(i,1) == 8888
                L(i,1) = 0.0544 + 0.7536*w(i,1) + 0.0297*w(i,1)^2;
                M(i,1) = 0.8678 -0.1785*w(i,1) + 0.1401*w(i,1)^2;
                N(i,1) = 2;
            end
        end
        alpha(i,1) = Tr(i,1)^(N(i,1)*(M(i,1)-1))*exp(L(i,1)*(1-Tr(i,1)^(M(i,1)*N(i,1))));
    case 'SOAVE'
        m = 0.37464 + 1.54226.*w - 0.26992.*w.^2;
        alpha = (1 + m.*(1 - sqrt(Tr))).^2;
end
% atractive parameter a(T) in bar . L² / mol²
ac = (0.45724.*((R.*Tc).^2)./Pc);
a = ac.*alpha;

% co-volume (b) in L/mol
b = 0.0778.*(R.*Tc)./Pc;

% Volume translation (c)
F = zeros(NC,1);
for i = 1:NC
    if c(i,1) == 8888
        F(i,1) = -0.014471 + 0.067498*w(i,1) - 0.084852*w(i,1)^2 + 0.067298*w(i,1)^3 - 0.017366*w(i,1)^4;
        c(i,1) = (R*Tc(i,1)/Pc(i,1))*F(i,1);
    end
end

% a/bRT
abRT = a./(b.*R.*T);

%==========================================================================
end