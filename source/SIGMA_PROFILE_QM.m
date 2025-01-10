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
% --> molecule = molecule's name [CHAR 1x1]
% --> method = COSMO or GCM (group contribution method) [CHAR 1x1]
%
% OUTPUTS
% --> p_sigma = charge density distribution vector [61x2]
%               - firts column = sigma-values in e/nm^2
%               - second column = p(sigma values) in nm^2
%
%==========================================================================

function p_sigma = SIGMA_PROFILE_QM(molecule)

filename = strcat('sigma-profiles\',molecule,'.dat');

if isfile(filename)
    sigma_profile = readtable(filename);
    p_sigma(:,1) = sigma_profile.Var1;
    p_sigma(:,2) = sigma_profile.Var2;
else
    p_sigma = zeros(61,2);
end

%--------------------------------------------------------------------------
end % function