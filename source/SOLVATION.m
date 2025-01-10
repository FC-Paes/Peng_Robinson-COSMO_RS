%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function calculates the Gibbs energy of solvation from the combined
% knowledge of solvent's dendity and solute's fugacity calculated by the 
% PR/COSMO-RS EoS
%==========================================================================
%
% INPUTS:
% -- T = temperature in K [REAL 1x1]
% -- P = pressure in bar [REAL 1x1]
% -- z = mixture composition [REAL NCx1]
% -- phase = phase to performo solvation free energy calculations
%            phase = 0 (most stable option)
%            phase = 1 (liquid)
%            phase = 2 (vapor)
% -- molecule_list = list with the COSMO name of each molecule [REAL NCx1]
% -- data = structure containing molecules data (critical properties, sigma, profiles, etc.)
% -- par = structure containing the parametrization of COSMO-RS
%
% OUTPUTS: 
% -- DG = Gibbs energy of solvation
%
% PS:  in the molecule_list, we provide:
% Firts row = solvent
% Second row = solute
%
%==========================================================================

function DG = SOLVATION(T,P,z,phase,molecule_list,par)

%--------------------------------------------------------------------------
% ln of fugacity coefficients and Volume of the liquid phase
[lnFUG,V] = THERMO(T,P,z,molecule_list,phase,par);

%--------------------------------------------------------------------------
% Solvation energy

% gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2; 

% liquid-phase density
ro = 1/V;

% fugacities
FUG = exp(lnFUG);

% solvation gibbs energy in kcal/mol
DG = (8.314*T.*log(P.*FUG/(R*T*ro)))./4180;

end