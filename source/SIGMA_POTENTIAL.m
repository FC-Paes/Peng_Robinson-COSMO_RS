%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function calculates the sigma potential through the  
% sucessive substitution method
%==========================================================================
%
% INPUTS:
% --> profile = sigma-profile of pure compounds or mixtures, in nm^2
% --> Tk = Temperature in K
% --> par = structure containing the list of COSMO-RS universal parameters 
%
% OUTPUT: 
% --> potential = sigma potential at Tk, in kcal/mol/nm^2
%
%==========================================================================
function potential = SIGMA_POTENTIAL(profile,Tk,par)

% nomalization of liquid-phase sigma-profile
sum_profile = sum(profile,1);
psigma = profile./sum_profile;

% inputs from the sigma-profile for the successive substitution method
N = size(psigma,1); % number of points for the discretized profile

% Universal parameters of COSMO-RS
R = par.univ_muliq.R;
aeff = par.univ_muliq.aeff;

% -- Total contact energy matrix (Et = Em + Ehb)
Et = par.univ_muliq.Et;

% Defining matrix At and Mu
aeffRT = aeff/(R*Tk);
At = exp(aeffRT * Et);
Mu = zeros(N,N);
for i = 1:N
    for j = 1:N
        Mu(i,j) = psigma(j)./At(i,j);
    end
end

% Loop for the successive substitution method
error = 1e3;
Zold = ones(N,1);
while error > 1e-3 % error in kcal/mol/nm^2
    Z = Mu * Zold;
    Z = 1./Z;
    error = max(abs(Z-Zold));
    Zold = (Z + Zold) / 2;
end
aeffRT_inv = (R * Tk) / aeff;
potential = aeffRT_inv * log(Z);
potential = 100.*potential./4.18; % [kcal/mol/nm^2]

end % function