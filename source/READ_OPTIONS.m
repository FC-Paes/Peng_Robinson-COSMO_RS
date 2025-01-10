%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% Reading calculation settings from '.\INPUTS.xlsx'
%==========================================================================
%
% INPUTS: 
%
% OUTPUT: 
% -- calc = simulation or optimization + simulation
% -- dataT = type of solvation data to be load (molecules, free radicals,
%            transition states, etc...)
% -- psgima_ipt = type of sigma profile to be used (QM-based or GC-based)
% -- level = level of theory of sigma profiles (only BP-TZVPD-FINE is
%            available in this implementation)
% -- alpha_function = T-dep alpha function to be used with the PR EoS
%            (Soave 1972 or Twu-91)
% -- EoS_par = method to calculate the PR EoS parameters (using Exp data or
%            GC-based methods)
%
%==========================================================================

function [calc,dataT,psgima_ipt,level,alpha_function,EoS_par] = READ_OPTIONS

options = readtable(".\INPUTS.xlsx",Sheet="Options");

calc = char(table2array(options(1,2)));
dataT = char(table2array(options(2,2)));
level = char(table2array(options(3,2)));
psgima_ipt.solvent = string(table2array(options(4,2)));
psgima_ipt.solute = string(table2array(options(5,2)));
EoS_par.solvent = string(table2array(options(6,2)));
EoS_par.solute = string(table2array(options(7,2)));
alpha_function = char(table2array(options(8 ,2)));

end