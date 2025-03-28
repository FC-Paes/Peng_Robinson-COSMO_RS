%==========================================================================
% function to load properties of molecules, such as:
% -- Tc, Pc, omega
% -- L,M,N for the twu-91 alpha function
% -- COSMO cavity data (Area, Volume, Edisp, sigma-profile,...)
% -- UNIFAC decomposition
%==========================================================================
%
% INPUTS
% -- ID = list of molecule - string [Nmol x 1]
% -- IDtype = type of information that was given in the ID vector (name, cas, or smiles) - string [1 x 1]
%
% OUTPUTS
% -- molecule_list = structure containing all the information of the molecules informed in ID - - stucture [Nmol x 1]
%
%==========================================================================
function molecule_list = DATA_MOL(ID,IDtype,data,methodPR,methodCOSMO,par)

%-----------------------------------
% PR EoS parameters

molecule_list.calc = methodPR;
switch methodPR
    case "GCM"
        switch IDtype
            case "name"
                data_molecule = data.gc(strcmp(data.gc.Name_simulis_,ID),:);
            case "cas"
                data_molecule = data.gc(strcmp(data.gc.CAS_RN,ID),:);
            case "smiles"
                data_molecule = data.gc(strcmp(data.gc.SMILES,ID),:);
        end

        % ID
        molecule_list.name = string(data_molecule.Name_simulis_(1));
        molecule_list.cas = string(data_molecule.CAS_RN(1));
        molecule_list.family = string(data_molecule.Family_simulis_(1));
        
        % functional groups decomposition
        molecule_list.groups = string(data_molecule.UNIFAC(1));

        % PR EoS parameters
        [Tc,Pc,w,ac,b,m] = CRITICAL_PROPERTIES_GC(molecule_list.groups,par);
        
        % add PR EoS parameters to the molecule list
        molecule_list.ac = ac;
        molecule_list.b = b;
        molecule_list.m = m;
        molecule_list.Tc = Tc;
        molecule_list.Pc = Pc;
        molecule_list.w = w;
        molecule_list.p_sigma = SIGMA_PROFILE_GC(molecule_list.groups,par);

    case "Exp data"
        switch IDtype
            case "name"
                data_molecule = data.exp(strcmp(data.exp.Name_simulis_,ID),:);
            case "cas"
                data_molecule = data.exp(strcmp(data.exp.CAS_RN,ID),:);
            case "smiles"
                data_molecule = data.exp(strcmp(data.exp.SMILES,ID),:);
        end
        % ID
        molecule_list.name = string(data_molecule.Name_simulis_(1));
        molecule_list.cas = string(data_molecule.CAS_RN(1));
        molecule_list.family = string(data_molecule.Family_simulis_(1));

        % EoS parameters
        molecule_list.Tc = data_molecule.Tc_K_(1); % [K]
        molecule_list.Pc = data_molecule.Pc_bar_(1); % [bar]
        molecule_list.w = data_molecule.w___(1); % [-]
        molecule_list.L = data_molecule.L___(1); % [-]
        molecule_list.M = data_molecule.M___(1); % [-]
        molecule_list.N = data_molecule.N___(1); % [-]
        molecule_list.c = data_molecule.c_cm3_mol_(1)/1000; % [L/mol]
end

%-----------------------------------
% COSMO-RS sigma-profiles

switch methodCOSMO
    case "GCM"
        molecule_list.p_sigma = SIGMA_PROFILE_GC(molecule_list.groups,par);

    case "COSMO"
        molecule_list.p_sigma = SIGMA_PROFILE_QM(molecule_list.name); % [e/nm2]
end

%-----------------------------------
% Defining Association Code using sigma-moments
[Mhb,Mi] = MOMENTS(molecule_list.p_sigma);
molecule_list.Mhb = Mhb;
molecule_list.Mi = Mi;

if Mhb(1) > 1e-3 && Mhb(2) > 1e-3
    molecule_list.AssociationCode = "SA";

elseif Mhb(1) > 1e-3 && Mhb(2) < 1e-3
    molecule_list.AssociationCode = "HA";

elseif Mhb(1) < 1e-3 && Mhb(2) > 1e-3
    molecule_list.AssociationCode = "HD";

else
    molecule_list.AssociationCode = "NA";

end

%-----------------------------------
end
