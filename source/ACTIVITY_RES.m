%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine

%==========================================================================
% This function calculates the  residual contribution to
% the activity coefficient of a compound X in a liquid phase
%==========================================================================
%
% INPUTS:
% -- Tk = temperature in K [REAL 1x1]
% -- molecule_list = list with the COSMO name of each molecule [REAL NCx1]
%                    - first column = names
%                    - second column = type 
%                    - third column = method for inputs (COSMO or GCM)
% -- properties = structure containing critical properties
% -- par = structure containing the parametrization of COSMO-RS
%
% OUTPUT: 
% --> lnGam_r = residual contribution to the activity coefficiente of the compound X in the mixture [REAL NCx1]
% --> muRES_mix = residual contribution to the chemical potential of the compound X in the mixture [REAL NCx1]
% --> muRES_ideal = residual contribution to the chemical potential of the compound X in itself [REAL NCx1]
%
%==========================================================================
function [lnGam_r,muRES_mix,muRES_ideal] = ACTIVITY_RES(z,Tk,molecule_list,par)

%--------------------
% number of compounds
NC = size(molecule_list,1);

%--------------------
% Sigma profile of pure compounds and mixture
profile = zeros(61,NC);
for i = 1:NC
    profile(:,i) = molecule_list{i,1}.p_sigma(:,2);
end
profile_mix = profile*z;

%--------------------
% SIGMA-POTENTIALS
% -- mixture
par.univ_muliq.Et = INTERACTION_MATRIX_TK(Tk,par.Em_298, par.Ehb_298, par);
potential_mix = SIGMA_POTENTIAL(profile_mix,Tk,par);
% -- pure compounds
potential_pure = 0.0*profile;
for i = 1:NC
    potential_pure(:,i) = SIGMA_POTENTIAL(profile(:,i),Tk,par);
end

%--------------------
% CHEMICAL POTENTIALS IN MIXTURE (molecule i in the mixture)
muRES_mix = potential_mix'*profile; % kcal/mol

%--------------------
% CHEMICAL POTENTIALS IDEAL SOLUTION (molecule i in pure i)
muRES_ideal = dot(potential_pure,profile); % kcal/mol

%--------------------
% ln(activity coefficient) - residual contribution
lnGam_r = ((muRES_mix - muRES_ideal)./(par.univ_muliq.R/4.184*Tk))';

end