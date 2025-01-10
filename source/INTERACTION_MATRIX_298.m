%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function generates the interaction matrix at 298 K
%==========================================================================
%
% INPUTS:
% --> par = list of COSMO-RS parameters 
%
% OUTPUT: 
% --> Em_298 = matrix containing the contact energy between segments (misfit interactions only)
% --> Ehb_298 = matrix containing the contact energy between segments  (hydrogen bond interactions only)
%
%==========================================================================

function [Em_298,Ehb_298] = INTERACTION_MATRIX_298(par)

% Universal parameters of COSMO-RS
sigma = [-0.03:0.001:0.03]'; %[e/A2]
N = length(sigma);
alpha = par.univ_muliq.alpha;
sigma_hb = par.univ_muliq.sigma_hb;
Chb = par.univ_muliq.Chb;

% Interaction matrix
Em_298 = zeros(N,N);
Ehb_298 = zeros(N,N);
for i = 1:N
    for j = 1:N
        % Interaction energies between sigma_i and sigma_j
        % -- Misfit (or electrostatic) energy
        Em_298(i,j) = Emisfit(sigma(i),sigma(j),alpha);
        % -- Hydrogen-bond energy
        Ehb_298(i,j) = Ehydbon(sigma(i),sigma(j),sigma_hb,Chb);
    end
end
%==========================================================================
% NESTED FUNCTIONS 
%==========================================================================
    function Emis = Emisfit(sigma_i,sigma_j,alpha)
        % Misfit (or electrostatic) energy
        Emis = (alpha/2)*(sigma_i + sigma_j)^2; % [J/mol/A2]
    end
%==========================================================================
    function Ehb = Ehydbon(sigma_i,sigma_j,sigma_hb,Chb)
        % -- Hydrogen-bond energy  
        % HB donnor
        sigma_don = min([sigma_i,sigma_j]);
        % HB acceptor
        sigma_acc = max([sigma_i,sigma_j]);
        % calculation of HB energy
        term1 = min([0,(sigma_don + sigma_hb)]);
        term2 = max([0,(sigma_acc - sigma_hb)]);
        Ehb = Chb*min([0,(term1*term2)]);
    end
%==========================================================================
end % function