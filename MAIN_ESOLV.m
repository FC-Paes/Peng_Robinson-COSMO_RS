%--------------------------------------------------------------------------
%
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
%
% This MATLAB script allows the user to perform calculations
% related to solvation data. The user can perform either simulation of
% solvation data or optimization + simulation. Note that the optimization 
% is done on to fine-tune the group parameters using solvation energies 
% as reference data.
%
% This script will:
% 1. Load solvation data based on molecule type (stable molecule, 
%    free radical, or transition state). Activation free energies of 
%    solvation can also be calculated for H-abstraction reactions.
% 2. Load input parameters, such as:
%    - Crical constants, acentric factors, Twu-91 alpha function parameters...
%    - COSMO-RS sigma-profiles
%    - Functional groups and their contribution parameters
%    - PR EoS and COSMO-RS universal constants
%    - etc...
% 3. Perform either a simulation or an optimization followed by simulation.
% 4. Generate plots (Parity Plot, Deviation Distribution Plot).
%
%------------
% PS: open the excel file 'INPUTS.xlxs' to define all the settings
%------------
%
%==========================================================================

clc,clear,close all
addpath(genpath("source\"))  

%-----------------------
% Load options

[calc,dataT,psigma_ipt,level,alpha_function,EoS_par] = READ_OPTIONS;

% ----------------------
% Load solvation data

switch dataT % data type to be considered...

    case 'Molecules' % stable solutes
        Esolv = readtable("solvation-data\ESOLV_DATA.xlsx","Sheet","Molecules");

    case 'Free radicals' % C/H/O radicals as solute
        Esolv = readtable("solvation-data\ESOLV_DATA.xlsx","Sheet","Free Radicals");

    case 'Transition states' % H-abstraction reactions
        Esolv = readtable("solvation-data\ESOLV_DATA.xlsx","Sheet","Transition states");

    case 'SFE of activation' % SFE = solvation free energy of activation related to H-abstraction reactions
        Esolv = readtable("solvation-data\ESOLV_DATA.xlsx","Sheet","SFE of activation");

    otherwise
        return

end

% ----------------------
% Load Parametrization

% -- PR/COSMO-RS universal parameters
par = PARAMETRIZATION_COSMO(level);
par.alpha_function = alpha_function;

% -- define the PR EoS parameters calculation method
% (using either experimental data or group contribution)
par.EoS_par = EoS_par; 

% -- define the sigma-profiles calculation methods
% (quantum-based profiles or calculated by group contribution)
par.psigma_ipt = psigma_ipt; 

% -- load functional groups and group contribution parameters
par.GCM = LOAD_GCM;

% ----------------------
% Load pure compound database
% (Critical constants, acentric factors, Twu 91 parameters, QM-based
% sigma-profiles, decomposition of compound into functional groups)

% Two tables are created to be used according to the input method chosen

% -- exp data for PR EoS parameters, along with quantum based sigma profile
data.exp = readtable("properties\Compounds_data.xlsm",Sheet='database_exp');
data.psigma_qm = readtable("sigma-profiles\p_sigma_res.xlsx",Sheet='QM',ReadVariableNames=true);

% -- molecules and group decomposition to apply the GC-based methods
data.gc = readtable("properties\Compounds_data.xlsm",Sheet='database_gc');

% ----------------------
% Calculations to be performed
% The calculation options are:
% -- simulation with the current parametrization of the GCM
% -- optimization of group contributions + simulation

switch calc

    case 'Simulation only'

        % simulation
        tic
        if dataT == "SFE of activation" 
            % solvation free energy of activation
            Esolv_res = SIMULATION_DDG(Esolv,data,par);

        else 
            % solvation free energy of molecules, free radicals or 
            % transition states
            Esolv_res = SIMULATION(Esolv,data,par);

        end
        toc

    case 'Optimization + simulation'

        % Here we load the database for optimization
        % and we chose a type of radical or TS to optimize

        Esolv = readtable("solvation-data\ESOLV_DATA.xlsx","Sheet","opt_data");
        molecules_types = ["all types";unique(Esolv.solute_type(:,1),'stable')];
        [indx,tf] = listdlg('ListString',molecules_types);
        type_for_opt = molecules_types(indx);

        if indx ~= 1
            Esolv = Esolv(strcmp(Esolv.solute_type,type_for_opt),:);
        end
        tic

        % optimization of associated group contribution parameters
        parNEW = OPTIMIZATION_GC(Esolv,data,par,type_for_opt);

        % simulation
        Esolv_res = SIMULATION(Esolv,data,parNEW);

        toc

    case 'None'
        return

end

% ----------------------
% Parity plots and deviations

Esolv_load = readtable("solvation-data\solv-energy-res.dat");
switch dataT
    case 'Molecules'
        PLOT_MOLECULES(Esolv_load)

    case 'Free radicals'
        PLOT_RADICALS(Esolv_load)

    case 'Transition states'
        PLOT_TS(Esolv_load)

    case 'SFE of activation'
        PLOT_ACTIVATION(Esolv_load)

end

% ----------------------
% JOB DONE =)
msgbox(["Simulation completed";"Check results in 'solvation-data/solv-energy-res.dat'"]);

