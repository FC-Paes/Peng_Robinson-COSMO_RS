%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
%
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% "SATURATION_PRESSURE_PURE" calculates the saturation pressure of a pure
% compound
%==========================================================================
%
% INPUTS
% -- T = temperature in K [REAL 1x1]
% -- Tc = critical temperature [REAL 1x1]
% -- Pc = Critical pressure [REAL 1x1] 
% -- w = acentric factor [REAL 1x1] 
% -- L, M, and N = Twu-91 parameters [REAL 1x1] 
% -- c = volume translation constant [REAL 1x1] 
% -- alpha_function = alpha function to be used (Twu-91 or Soave 1972)
%
% OUTPUTS
% -- Psat = saturation pressure in bar [REAL 1 x 1]
% -- lnFUGliq = ln(fugacity coefficient) of the liquid phase [REAL 1 x 1]
% -- lnFUGvap = ln(fugacity coefficient) of the vapor phase  [REAL 1 x 1]
%
%==========================================================================

function [Psat,lnFUGliq,lnFUGvap] = SATURATION_PRESSURE_PURE(T,Tc,Pc,w,L,M,N,c,alpha_function)

%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2;

% tolerance (stop criteria)
Tol = 1e-7;

%--------------------------------------------------------------------------
% Initialization for Psat using Vc
Vc = R.*Tc./Pc;
Tr = T./Tc;

% EoS parameters
[a,b,c,abRT] = PARAMETERS_EoS(T,Tc,Pc,w,L,M,N,c,alpha_function);

% Correcting P if P < 0 for Vm = Vc
P = pPR(a,b,c,T,Tr*Vc,R);
if P < 0
    k = 1;
    while P<=0
        Vm = Vc*(2^k);
        P = pPR(a,b,c,T,Vm,R);
        k = k + 1;
        if k > 1000
            P = 1;
            break
        end
    end
end

% last check before the main loop
% -- somtimes when P is negative for Vc, and we deduce Vm until P is
% becomes positives gives non-sense results, like P > 1e5 bar
% -- This come to happen to substances with high Tc and Pc

if P > 1e3
    P = 1; % New initialization
end

%--------------------------------------------------------------------------
% Full Calculation of Psat by Newton's method

% Main loop

Dev_psi = 1e3;
Dev_Psat = 1e3;
iter = 0;
while Dev_psi  > Tol || Dev_Psat > Tol
    iter = iter + 1;

    % calculation of molar volumes (liq and vap)
    [Zliq,Vliq] = MOLAR_VOLUME(T,P,b,abRT,1);
    [Zvap,Vvap] = MOLAR_VOLUME(T,P,b,abRT,2);

    % calculation of fugacities and Volumes at T and P (step k)
    lnFUGliq = fugPR(a,b,T,P,Zliq,R);
    lnFUGvap = fugPR(a,b,T,P,Zvap,R);

    % calculation of psi = lnFUGvap - lnFUGliq
    psi = lnFUGvap - lnFUGliq;
    P_new = P - R.*T.*(lnFUGvap - lnFUGliq)./(Vvap - Vliq);
    if P_new < 0
        P_new = P*exp(-(R.*T./P).*(lnFUGvap - lnFUGliq)./(Vvap - Vliq));
    end

    % chack for non-real results
    if isreal(P_new) == 0
        P_new = real(P_new); % bar
    end

    % stop criteria
    Dev_psi = max(abs(psi));
    Dev_Psat = abs(P_new - P)./P;

    % update value
    P = P_new;

end

Psat = P;

%==========================================================================

    function P = pPR(a_PR,b_PR,c_PR,Tk,Vm,R)
    % !  Pressure calculation by PR EoS
    % !  Inputs : a_PR, value of T-dep coefficient in bar.L^2/mol^2[REAL]
    % !           b_PR, covolume of the mixture in L/mol [REAL]
    % !           tk, temperature in K
    % !  Output : pPR in Bar

    P = R*Tk/(Vm + c_PR - b_PR) - a_PR/((Vm + c_PR)*(Vm + c_PR + b_PR) + b_PR*(Vm + c_PR - b_PR));
    end

%==========================================================================

    function phi = fugPR(a_PR,b_PR,Tk,Pbar,Z,R)
    % !  Fugacity coefficient calculation by PR EoS
    % !  Inputs : a_PR, value of T-dep coefficient in bar.L^2/mol^2[REAL]
    % !           b_PR, covolume of the mixture in L/mol [REAL]
    % !           tk, temperature in K
    % !  Output : pPR in Bar
    %-------------

    % EoS parameters
    A = a_PR*Pbar/(R*Tk)^2;
    B = b_PR*Pbar/(R*Tk);

    % ln of the fugacity coefficient (pure compound)
    phi = (Z - 1 - log(Z-B) - A/(2*B*sqrt(2))*log((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B)));

    end
%==========================================================================
end % function
