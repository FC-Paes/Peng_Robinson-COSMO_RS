%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
%
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function calculates the fugacity coefficient of each compound in the mixture,
% as well as the molar volume of the bulk phase.
%==========================================================================
%
% INPUTS:
% -- T = temperature in K [REAL 1x1]
% -- P = pressure in bar [REAL 1x1]
% -- z = mixture composition [REAL NCx1]
% -- molecule_list = list with the COSMO name of each molecule [REAL NCx1]
% -- phase = 1 is liquid, 2 is vapor
% -- par = structure containing the parametrization of the PR/COSMO-RS EoS
%
% OUTPUTS:
% -- lnFUG = ln of the fugacity coefficients [REAL NCx1]
% -- V = mixture molar volume in L/mol [REAL 1x1]
%
%==========================================================================

function [lnFUG,V] = THERMO(T,P,z,molecule_list,phase,par)

%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2; 

%--------------------------------------------------------------------------
% Activity coefficients calculation with COSMO-RS

% number of compounds
NC = size(molecule_list,1);

% residual contribution from COSMO-RS
lnGam_r = ACTIVITY_RES(z,T,molecule_list,par);

% combinatorial part
lnGam_c = 0; % zero for a HV-based mixing rule (Pref = +inf)

% all contributions
lnGam = lnGam_r + lnGam_c;

%--------------------------------------------------------------------------
% EoS parameters

b = zeros(NC,1);
c = zeros(NC,1);
abRT = zeros(NC,1);
for i = 1:NC

    if molecule_list{i,1}.calc == "Exp data"

        % Inputs required to calculate the PR EoS parameters
        Tc = molecule_list{i,1}.Tc;
        Pc = molecule_list{i,1}.Pc;
        w = molecule_list{i,1}.w;
        L = molecule_list{i,1}.L;
        M = molecule_list{i,1}.M;
        N = molecule_list{i,1}.N;
        TT = molecule_list{i,1}.c;

        % EoS parameters calculation
        if strcmp(molecule_list{i}.name,"h2")==1 || strcmp(molecule_list{i}.name,"he")==1 % quantum fluids
            [~,b(i,1),c(i,1),abRT(i,1)] = PARAMETERS_EoS(T,Tc,Pc,w,L,M,N,TT,'SOAVE');

        else % other fluids we can use either Twu-91 or Soave 1972 as alpha-function
            [~,b(i,1),c(i,1),abRT(i,1)] = PARAMETERS_EoS(T,Tc,Pc,w,L,M,N,TT,par.alpha_function);
        end

    elseif molecule_list{i,1}.calc == "GCM"

        % The parameters of the PR EoS are calculated based on group
        % contribution results
        ac = molecule_list{i,1}.ac;
        b(i,1) = molecule_list{i,1}.b;
        m = molecule_list{i,1}.m;
        Tr = T/molecule_list{i,1}.Tc;
        alpha = (1 + m*(1 - sqrt(Tr)))^2;
        abRT(i,1) = (ac*alpha)/(b(i,1)*R*T);
        c(i,1) = 0.0;

    end
end

%--------------------------------------------------------------------------
% Mixing rules to calculate the mixture parameters
[b_mix,c_mix,abRT_mix,sum_zb,der_abRT_mix] = MIX_RULE_EOS(z,b,c,abRT,lnGam,par);

%--------------------------------------------------------------------------
% Here we solve the PR EoS for V
% phase = 1 (liquid)
% phase = 2 (vapor)
[~,V] = MOLAR_VOLUME(T,P,b_mix,abRT_mix,phase);

%--------------------------------------------------------------------------
% Now we calculate ln(Fugacity coefficients) of all compounds in the
% mixture
lnFUG = FUGACITY(T,P,V,b_mix,sum_zb,der_abRT_mix,c);

%--------------------------------------------------------------------------
% Finaly we correct the molar volume using the volume translation constant
V = V - c_mix;

end