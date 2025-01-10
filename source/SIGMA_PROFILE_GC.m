%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
%
% This function searches the sigma-profile of a given molecule in the
% database
%
%==========================================================================
% INPUTS
% -- molecule = molecule's name [CHAR 1x1]
%
% OUTPUTS
% -- p_sigma = charge density distribution vector [61x2]
%               - firts column = sigma-values in e/nm^2
%               - second column = p(sigma values) in nm^2
%
%==========================================================================
function p_sigma = SIGMA_PROFILE_GC(groups_string,par)

% Here we perform the complete calculation of the sigma profile (by
% group contribution), that can be directly used in the calculation
% of activity coefficients
% occurances
groups_UNIFAC = par.GCM.groups_UNIFAC;
occ = DECOMPOSITION_VECTOR(groups_string,groups_UNIFAC);
% parameters
b_sp = par.GCM.sigma_profile;
% sigma-profile
p_sigma(:,1) = par.GCM.sigma_values;
p_sigma_raw = occ*b_sp;
p_sigma_raw(p_sigma_raw < 0) = 0;
p_sigma(:,2) = p_sigma_raw';

end
