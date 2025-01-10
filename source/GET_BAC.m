%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% Function used to obtain the binary association (BAC) code of a given,
% based on the association nature of its compounds, which can be:
% NA - non-associating
% HA - hydrogen-bonding donor
% HD - hydrogen-bonding acceptor
% SA - self-associating
%==========================================================================
%
% INPUTS: 
% -- AC_1 = association nature of compound 1
% -- AC_2 = association nature of compound 1
%
% OUTPUT: 
% -- BAC = binary association code
% -- Type = Type of BAC:
%                      1 - association does not take place (BAC 1-4)
%                      2 - only self association takes place (BAC 5)
%                      3 - only cross association takes place (BAC 6)
%                      4 - cross and self association take place (BAC 7-9)
%
%==========================================================================

function [BAC,Type] = GET_BAC(AC_1,AC_2)

CODE = append(AC_1,'-',AC_2);
    % Binary association code
    switch CODE
        case 'NA-NA' % category 1
            BAC = 1;
            Type = 1;
        case {'NA-HA','HA-NA'} % category 2
            BAC = 2;
            Type = 1;
        case {'HD-NA','NA-HD'} % category 3
            BAC = 3;
            Type = 1;
        case {'HA-HA','HD-HD'} % category 4
            BAC = 4;
            Type = 1;
        case {'SA-NA','NA-SA'} % category 5
            BAC = 5;
            Type = 2;
        case {'HD-HA','HA-HD'} % category 6
            BAC = 6;
            Type = 3;
        case {'SA-HD','HD-SA'} % category 7
            BAC = 7;
            Type = 4;
        case {'SA-HA','HA-SA'} % category 8
            BAC = 8;
            Type = 4;
        case {'SA-SA'} % category 9
            BAC = 9;
            Type = 4;
        otherwise % nonassessed cactegory
            BAC = 0;
            Type = 0;
    end