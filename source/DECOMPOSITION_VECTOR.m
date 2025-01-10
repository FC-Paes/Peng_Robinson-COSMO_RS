%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function gives the occurance vector of UNIFAC groups, based on a 
% string containing the groups of the molecule. 
% Ex: Propane:  [CH3]: 2, [CH2]: 1
%     Propanol: [CH3]: 1, [CH2]: 1, [OH]: 1
%==========================================================================
%
% INPUTS:
% -- groups_string = string with UNIFAC groups of a given molecule
% -- groups_UNIFAC = list of all groups
%
% OUTPUT: 
% --> occ = occurance vector
%
%==========================================================================
function occ = DECOMPOSITION_VECTOR(groups_string,groups_UNIFAC)

groups_string = strcat(groups_string,',');
N = size(groups_UNIFAC,1);
occ = zeros(1,N);
for i = 1:N
    group = (groups_UNIFAC(i));
    if contains(groups_string,char(group))
        newStr1 = extractAfter(groups_string,strcat(char(group),':'));
        newStr2 = extractBefore(newStr1,',');
        occ(i) = str2double(newStr2);
    else
        occ(i) = 0;
    end
end

end