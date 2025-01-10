%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function sets the GC method:
% It loads an excel file containing the parameters of all groups
% for the calculation of :
% -- sigma-profiles
% -- PR EoS parameters (ac, b, and a0 - see main text)
%
%==========================================================================
%
% INPUTS:
%
% OUTPUT: 
% --> GCM = structure contaning a decomposition list, as well as the
%           parameters of each group to predict cavity information
%
%==========================================================================
function GCM = LOAD_GCM

groups_UNIFAC = readtable("properties\Compounds_data.xlsm",Sheet="groups_UNIFAC");
GCM.groups_UNIFAC = table2array(groups_UNIFAC(:,2));
GCM.sigma_profile = readtable("properties\Compounds_data.xlsm",Sheet="psigma_GC",ReadVariableNames=true);
GCM.sigma_profile = table2array(GCM.sigma_profile(:,3:end));
GCM.sigma_values = [-3:0.1:3]';
EoS = readtable("properties\Compounds_data.xlsm",Sheet="EoS_GC",ReadVariableNames=true);
GCM.b = EoS.b;
GCM.ac = EoS.ac;
GCM.a0 = EoS.a0;

end % functions