%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function calculates the sigma-moments from sigma profiles
% This information is used to define the associative nature of molecules
% NA = Non-associating
% HA = Hydrogen-bonding acceptor
% HD = Hydrogen-bonding donor
% SA = Self-associating
%==========================================================================
%
% INPUTS:
% --> p_sigma = sigma-profile - 1st column = sigma values in e/nm2 [REAL 61 x 1]
%                             - 2nd column = probability values times cavity surface in nm2 [REAL 61 x 1]
% --> par = list of COSMO-RS parameters 
%
% OUTPUT: 
% --> Mhb = hydrogen bond sigma moments - 1st pos = HB-donor 
%                                       - 2nd pos = HB-acceptor
% --> Mi = sigma-moments from 0th to 7th order
%
%==========================================================================

function [Mhb,Mi] = MOMENTS(p_sigma)

% function to calculate the sigma-moments for a given sigma profile
% input: p_sigma = sigma profile [61x2]
%                  Columns 1 = sigma values [e/Nm²]
%                  Columns 1 = probability values x Total surface [Nm²]
% outputs: Mi = moments of the given sigma profile distribution
s = p_sigma(:,1);
p = p_sigma(:,2);

% sigma-moments initialization
Mi = zeros(1,7);
Mhb = zeros(1,2);


% sigma-moments calculation (M0,M1,...,M6)
for i = 1:7
    m = i-1;
    fm = s.^m;
    Mi(i) = dot(p,fm);
end

% sigma-moments calculation (Macc, Mdon)
shb = 1.1; %[e/Nm²] = cut-off value for HB
facc= zeros(size(s,1),1);
fdon= zeros(size(s,1),1);
for i = 1:size(s,1)
    % (HB acc)
    if s(i) > shb
        facc(i,1) = s(i) - shb;
    end
    % (HB don)
    if -s(i) > shb
        fdon(i,1) = -s(i) - shb;
    end
end
Mhb(1) = dot(p,facc); % (HB acc)
Mhb(2) = dot(p,fdon); % (HB acc)