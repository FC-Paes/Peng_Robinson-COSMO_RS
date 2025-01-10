%==========================================================================
% "CRITICAL_PROPERTIES" get the critical and COSMO properties of stable molecules
%==========================================================================
%
% INPUTS
% -- molecule_list = list with the COSMO name of each molecule [REAL NCx1]
%                    - first column = names
%                    - second column = type (stable or radical)
%                    - third column = method for inputs (COSMO or GCM)
% -- properties = structure containing critical properties
%
% OUTPUTS
% -- Tc = vector of critical temperatures in K [REAL NC x 1]
% -- Pc = vector of critical pressures in bar [REAL NC x 1]
% -- w = vector of acentric vectors [REAL NC x 1]
% -- Area = vector of cavity surfaces in nm^2 [REAL NC x 1]
% -- Vol = vector of cavity volumes in nm^3 [REAL NC x 1]
% -- Eint = vector of pure compound interaction energy [REAL NC x 1]
%
%==========================================================================

function [Tc,Pc,w,ac,b,m] = CRITICAL_PROPERTIES_GC(groups_string,par)

%  gas constant [=] bar.L/(mol K)
R = 8.31446261815324e-2;

%---------------------------
% Model constants
omega_a = 0.45724;
omega_b = 0.0778;
expb = 4/5;
expac = 2/3;
expaT0 = 0.48;
ctes_b = par.GCM.b;
ctes_ac = par.GCM.ac;
ctes_aT0 = par.GCM.a0;

OmhR_a = omega_a*R^2;
OmhR_b = omega_b*R;

%---------------------------
% Groups availables
groups_UNIFAC = par.GCM.groups_UNIFAC;

%---------------------------
% decomposition vector
occ = DECOMPOSITION_VECTOR(groups_string,groups_UNIFAC);

%---------------------------
% Parameters EoS
b = OmhR_b*(([occ]*ctes_b).^(1/expb));
ac = OmhR_a*(([occ]*ctes_ac).^(1/expac));
aT0 = OmhR_a*(([occ]*ctes_aT0).^(1/expaT0));

%---------------------------
% critical temperature, Tc
Tc = (omega_b/omega_a*1/R).*(ac./b);

%---------------------------
% Critical pressure, Pc
Pc = (omega_b^2/omega_a).*(ac./(b.^2));

%---------------------------
% Parameter m for the Soave alpha function
m = (((aT0))./(ac)).^0.5 - 1;

%---------------------------
% acentric factor

% -- Coefficients of the quadratic equation
a1 = -0.26992;
b1 = 1.5422;
c1 = 0.37464 - m; % Subtract m to bring equation to standard form

% Compute discriminant
D = b1^2 - 4*a1*c1;

% Check if the discriminant is non-negative
if D >= 0

    % Compute the two solutions
    w1 = (-b1 + sqrt(D)) / (2*a1);
    w2 = (-b1 - sqrt(D)) / (2*a1);

    if w1 > -0.5 && w1 < 3
    % if only w1 is within a reasonable range for acentric factors, 
    % we chose w1
        w = w1;
    elseif w2 > -0.5 && w2 < 3
    % if only w2 is within a reasonable range for acentric factors, 
    % we chose w2
        w = w2;
    elseif w1 > -0.5 && w1 < 3 && w2 > -0.5 && w2 < 3.0
    % if both w1 and w2 are within a reasonable range for acentric factors, 
    % we chose the one closest to zero (tipically acentric factors are
    % small values
        w = min([abs(w1),abs(w2)]);
    else
    % otherwise, we assume w = 0 to avoid abnormal results
        w = 0;
    end
    
else
    % Handle the case of no real solutions
    w = 0;
end

%==========================================================================
end