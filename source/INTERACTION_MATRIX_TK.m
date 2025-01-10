%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function corrects the interaction matrix from 298 K to other
% temperatures
%==========================================================================
%
% INPUTS:
% --> Tk = matrix containing the misfits energy between segments at 298 K
% --> Em_298 = matrix containing the misfits energy between segments at 298 K
% --> Ehb_298 = matrix containing the HB energy between segments at 298 K
%
% OUTPUT: 
% --> Em_tk = matrix containing the misfits energy between segments at Tk
% --> Ehb_tk = matrix containing the HB energy between segments at Tk
%
%==========================================================================

function [Et,Em,Ehb] = INTERACTION_MATRIX_TK(Tk,Em_298, Ehb_298, par)

% -- Hydrogen-bond energy  
% Temperature dependence of the Hydrogen-bond constant 
ChbT1 = par.univ_muliq.ChbT1;
ChbT2 = par.univ_muliq.ChbT2;
termT = Tk*log(1 + exp(ChbT1*1e3/(8.314*Tk)) * abs(ChbT2));
term298 = 298.15*log(1 + exp((ChbT1*1e3)/(8.314*298.15)) * abs(ChbT2));
fhbT = termT/term298;
Ehb = fhbT * Ehb_298;


% Temperature dependence of the Hydrogen-bond constant 
% -- Misfit energy  
fmT = 1.0; % no T-function added in this version of COSMO-RS
Em = fmT * Em_298;

% Final interaction matrix
Et = Em + Ehb;

%==========================================================================
end % function













