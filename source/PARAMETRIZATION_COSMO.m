%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function loads the universal parameters of COSMO-RS
% as well as the contants of the mixing rule
%==========================================================================

function par = PARAMETRIZATION_COSMO(level)

% Constants to calculate the chemical potential in liquid phase with
% COSMO-RS

ctesLIQ = readtable(".\INPUTS.xlsx",Sheet="Parametrization_COSMO");
switch level
    case "BP-TZVP"
        par.univ_muliq.aeff = ctesLIQ.aeff_A_2_(1,1); % A^2
        par.univ_muliq.alpha = ctesLIQ.alpha_kJ_mol_A_2_(1,1);% [kJ/mol/A^2]
        par.univ_muliq.Chb = ctesLIQ.Chb_kJ_mol_A_2_(1,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT1 = ctesLIQ.ChbT1___(1,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT2 = ctesLIQ.ChbT2___(1,1); % [kJ/mol/A^2]
        par.univ_muliq.sigma_hb = ctesLIQ.sigma_hb_e_A_2_(1,1); % e/A^2
        par.univ_muliq.Cvdw = ctesLIQ.Cvdw___(1,1); % [kJ/mol/A^2]
        par.univ_muliq.CvdwT = ctesLIQ.CvdwT___(1,1); % [kJ/mol/A^2]
        par.univ_muliq.R = 8.31446261815324/1000; % kJ/mol/K
    case "BP-TZVPD-FINE"
        par.univ_muliq.aeff = ctesLIQ.aeff_A_2_(2,1); % A^2
        par.univ_muliq.alpha = ctesLIQ.alpha_kJ_mol_A_2_(2,1);% [kJ/mol/A^2]
        par.univ_muliq.Chb = ctesLIQ.Chb_kJ_mol_A_2_(2,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT1 = ctesLIQ.ChbT1___(2,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT2 = ctesLIQ.ChbT2___(2,1); % [kJ/mol/A^2]
        par.univ_muliq.sigma_hb = ctesLIQ.sigma_hb_e_A_2_(2,1); % e/A^2
        par.univ_muliq.Cvdw = ctesLIQ.Cvdw___(2,1); % [kJ/mol/A^2]
        par.univ_muliq.CvdwT = ctesLIQ.CvdwT___(2,1); % [kJ/mol/A^2]
        par.univ_muliq.R = 8.31446261815324/1000; % kJ/mol/K
    case "DMOL3"
        par.univ_muliq.aeff = ctesLIQ.aeff_A_2_(3,1); % A^2
        par.univ_muliq.alpha = ctesLIQ.alpha_kJ_mol_A_2_(3,1);% [kJ/mol/A^2]
        par.univ_muliq.Chb = ctesLIQ.Chb_kJ_mol_A_2_(3,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT1 = ctesLIQ.ChbT1___(3,1); % [kJ/mol/A^2]
        par.univ_muliq.ChbT2 = ctesLIQ.ChbT2___(3,1); % [kJ/mol/A^2]
        par.univ_muliq.sigma_hb = ctesLIQ.sigma_hb_e_A_2_(3,1); % e/A^2
        par.univ_muliq.Cvdw = ctesLIQ.Cvdw___(3,1); % [kJ/mol/A^2]
        par.univ_muliq.CvdwT = ctesLIQ.CvdwT___(3,1); % [kJ/mol/A^2]
        par.univ_muliq.R = 8.31446261815324/1000; % kJ/mol/K
end

% Mixing rule constants for the EoS/gE approach
ctesMIX = readtable(".\INPUTS.xlsx",Sheet="Parametrization_Mix_Rule");
par.univ_mixrule.q1 = ctesMIX.q1___(1,1); 
par.univ_mixrule.s = ctesMIX.s___(1,1);

end
%==========================================================================