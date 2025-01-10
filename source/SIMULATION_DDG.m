%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
%
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
% This function performs solvation free energy calculations from data
% informed in the input file "\solvation-data\ESOLV_DATA.xlsx"
%==========================================================================
%
% INPUTS:
% -- Esolv = table containing information about solvent, solute,
%            temperature, pressure, and solvation free energies
% -- data = structure containing pure compound data, such as: critical
%           constants, sigma profiles, smiles, group decomposition, etc.
% -- par = structure containing the parametrization of the PR/COSMO-RS EoS
%
% OUTPUTS:
% -- Esolv_res = same as Esolv but with solvation prediction results
% (additional column)
%
%==========================================================================

function Esolv_res = SIMULATION_DDG(Esolv,data,par)

%--------------------------------------------------------------------------
% file head (answer file)
fileID = fopen('solvation-data\solv-energy-res.dat','w');
fprintf(fileID,'%150s %50s %100s %150s %100s %150s %100s %150s %100s %30s %30s %30s %30s %30s %30s %30s %30s %30s %30s\n',...
    'Solvent','CAS_solvent','Family_solvent',...
    'Reactant1','Family_R1',...
    'Reactant2','Family_R2',...
    'Transition-State','Family_TS',...
    'T[K]','P[bar]',...
    'DGsolv_exp_1[kcal/mol]','DGsolv_calc_1[kcal/mol]',...
    'DGsolv_exp_2[kcal/mol]','DGsolv_calc_2[kcal/mol]',...
    'DGsolv_exp_TS[kcal/mol]','DGsolv_calc_TS[kcal/mol]',...
    'DDGsolv_exp[kcal/mol]','DDGsolv_calc[kcal/mol]');

% Number of exp points to be simulated
Npoints = size(Esolv,1);

% wait bar to follow the progress of the calculation
msg = strcat('calculation progress:',num2str(0),' out of ',num2str(Npoints));
f = waitbar(0,msg,'Name','Calculation of solvation free energies');

%--------------------------------------------------------------------------
% Main loop 
Esolv_res = Esolv;

% Initialization of the interaction matrix at 298 K
[par.Em_298,par.Ehb_298] = INTERACTION_MATRIX_298(par);

for i = 1:Npoints
    %-------------------------
    % update of the wait bar
    msg = strcat('calculation progress:',num2str(i),'-->',num2str(Npoints));
    waitbar(i/Npoints,f,sprintf(msg))
    %-------------------------
    % Point information (solvent, solute, temperature, and pressure)
    solvent = string(Esolv.solvent(i));
    solute1 = string(Esolv.solute1(i));
    solute2 = string(Esolv.solute2(i));
    solute3 = string(Esolv.solute3(i));
    Tk = Esolv.T_K_(i);
    isHead = ismember('P_bar_', Esolv.Properties.VariableNames);
    if isHead
        Pbar = Esolv.P_bar_(i);
    else
        Pbar = 0; % in this case we are going to use the saturation pressure
    end
    DGsolv_exp_1 = Esolv.DGsolv1_kcal_mol_(i);
    DGsolv_exp_2 = Esolv.DGsolv2_kcal_mol_(i);
    DGsolv_exp_TS = Esolv.DGsolvTS_kcal_mol_(i);
    %-------------------------
    % check if both solvent and solute are present in the database
    % (for both critical properties and sigma-profiles)

    % solvent
    % --> solvent
    stat_solvent = 0;

    %   - check PR inputs
    if par.EoS_par.solvent == "Exp data"
        % check if PR EoS inputs are available in the experimental database
        stat_solvent_crit = ismember(solvent,data.exp.Name_simulis_);
    else
        % Ohterwise we check the group decomposition
        stat_solvent_crit = ismember(solvent,data.gc.Name_simulis_);
    end

    %   - check sigma profiles
    if par.EoS_par.solvent == "Exp data"
        % check if sigma-profiles are available in the quantum-based database
        stat_solvent_sp = 0;
        filename = strcat('sigma-profiles\',solvent,'.dat');
        if isfile(filename)
            stat_solvent_sp = 1;
        end
    else
        % Ohterwise we check the group decomposition
        stat_solvent_sp = ismember(solvent,data.gc.Name_simulis_);
    end

    % Can the solvent be used?
    % -- stat_solvent = 0 (No)
    % -- stat_solvent = 1 (yes)
    if stat_solvent_sp == 1 && stat_solvent_crit == 1
        stat_solvent = 1;
    end

    %===========
    % --> solute 1 (reactant 1 = stable molecule)
    stat_solute1 = 0;

    %   - check PR inputs
    if par.EoS_par.solute == "Exp data"
        % check if PR EoS inputs are available in the experimental database
        stat_solute1_crit = ismember(solute1,data.exp.Name_simulis_);
    else
        % Ohterwise we check the group decomposition
        stat_solute1_crit = ismember(solute1,data.gc.Name_simulis_);
    end

    %   - check sigma profiles
    if par.EoS_par.solute == "Exp data"
        % check if sigma-profiles are available in the quantum-based database
        stat_solute1_sp = 0;
        filename = strcat('sigma-profiles\',solute1,'.dat');
        if isfile(filename)
            stat_solute1_sp = 1;
        end
    else
        % Ohterwise we check the group decomposition
        stat_solute1_sp = ismember(solute1,data.gc.Name_simulis_);
    end

    % Can the solute be used?
    % -- stat_solute = 0 (No)
    % -- stat_solute = 1 (yes)
    if stat_solute1_sp == 1 && stat_solute1_crit == 1
        stat_solute1 = 1;
    end

    %===========
    % --> solute 2 (reactant 2 = free radical)
    stat_solute2 = 0;

    %   - check PR inputs
    if par.EoS_par.solute == "Exp data"
        % check if PR EoS inputs are available in the experimental database
        stat_solute2_crit = ismember(solute2,data.exp.Name_simulis_);
    else
        % Ohterwise we check the group decomposition
        stat_solute2_crit = ismember(solute2,data.gc.Name_simulis_);
    end

    %   - check sigma profiles
    if par.EoS_par.solute == "Exp data"
        % check if sigma-profiles are available in the quantum-based database
        stat_solute2_sp = 0;
        filename = strcat('sigma-profiles\',solute2,'.dat');
        if isfile(filename)
            stat_solute2_sp = 1;
        end
    else
        % Ohterwise we check the group decomposition
        stat_solute2_sp = ismember(solute2,data.gc.Name_simulis_);
    end

    % Can the solute be used?
    % -- stat_solute = 0 (No)
    % -- stat_solute = 1 (yes)
    if stat_solute2_sp == 1 && stat_solute2_crit == 1
        stat_solute2 = 1;
    end

    %===========
    % --> solute 3 (transition state)
    stat_solute3 = 0;

    %   - check PR inputs
    if par.EoS_par.solute == "Exp data"
        % check if PR EoS inputs are available in the experimental database
        stat_solute3_crit = ismember(solute3,data.exp.Name_simulis_);
    else
        % Ohterwise we check the group decomposition
        stat_solute3_crit = ismember(solute3,data.gc.Name_simulis_);
    end

    %   - check sigma profiles
    if par.EoS_par.solute == "Exp data"
        % check if sigma-profiles are available in the quantum-based database
        stat_solute3_sp = 0;
        filename = strcat('sigma-profiles\',solute3,'.dat');
        if isfile(filename)
            stat_solute3_sp = 1;
        end
    else
        % Ohterwise we check the group decomposition
        stat_solute3_sp = ismember(solute3,data.gc.Name_simulis_);
    end

    % Can the solute be used?
    % -- stat_solute = 0 (No)
    % -- stat_solute = 1 (yes)
    if stat_solute3_sp == 1 && stat_solute3_crit == 1
        stat_solute3 = 1;
    end

    % main loop
    if stat_solvent == 1 && stat_solute1 == 1 && stat_solute2 == 1 && stat_solute3 == 1
        % at this point the calculation can be performed
        %--------------------------
        IDtype = "name";

        % solvent
        molecule_list{1,1} = DATA_MOL(solvent,IDtype,data,par.EoS_par.solvent,par.psigma_ipt.solvent,par);

        % solute_1
        molecule_list{2,1} = DATA_MOL(solute1,IDtype,data,par.EoS_par.solute,par.psigma_ipt.solute,par);

        % solute_2
        molecule_list{3,1} = DATA_MOL(solute2,IDtype,data,par.EoS_par.solute,par.psigma_ipt.solute,par);

        % solute_3
        molecule_list{4,1} = DATA_MOL(solute3,IDtype,data,par.EoS_par.solute,par.psigma_ipt.solute,par);

        %--------------------------
        % -- composition
        % (1) = solvent
        % (2 to 4) = solutes
        z=[1;0;0;0];

        %--------------------------
        %--------------------------
        % Pressure check
        % for VERY low temperatures the Psat calculation bugs

        if Tk < 100
            % in such a case, we check if the pressure is given
            if Pbar > 0
                Pcalc = Pbar+0.2;
            else
                % otherwise we assume that the pressure is 1 bar
                Pcalc = 1; % bar
            end
        else

            % for temperatures above 100 K, we use the function to
            % calculate Psat of pure compounds

            if par.EoS_par.solvent == "Exp data"
                Tc_solv = molecule_list{1,1}.Tc;
                Pc_solv = molecule_list{1,1}.Pc;
                w_solv = molecule_list{1,1}.w;
                L_solv = molecule_list{1,1}.L;
                M_solv = molecule_list{1,1}.M;
                N_solv = molecule_list{1,1}.N;
                c_solv = molecule_list{1,1}.c;
                Psat = SATURATION_PRESSURE_PURE(Tk,Tc_solv,Pc_solv,w_solv,...
                    L_solv,M_solv,N_solv,c_solv,par.alpha_function);
            else

                groups_solv = molecule_list{1,1}.groups;
                [Tc_solv,Pc_solv,w_solv,~,~,~] = CRITICAL_PROPERTIES_GC(groups_solv,par);
                L_solv = 8888;
                M_solv = 8888;
                N_solv = 8888;
                c_solv = 8888;
                % note that 8888 is an error code, to point out that L, M, M, and c
                % are unknown, so generalized correlations can be used to
                % predict them based on Tc, Pc, and w
                Psat = SATURATION_PRESSURE_PURE(Tk,Tc_solv,Pc_solv,w_solv,...
                    L_solv,M_solv,N_solv,c_solv,par.alpha_function);

            end

            % Final pressure to be used in calculations
            % -- Psat or the pressure informed in the input file of
            % solvation data

            Pcalc = max(Psat,Pbar)+1e-6;

        end

        %--------------------------
        % calculation of solvation free energies
        % calc = 0 (GDsolv only)
        % calc = 1 (DGsolv + DSsolv + DHsolv)

        phase = 1; % liquid
        [DGsolv_calc] = SOLVATION(Tk,Pcalc,z,phase,molecule_list,par);
        DDG_exp = DGsolv_exp_TS - DGsolv_exp_2 - DGsolv_exp_1;
        DDG_calc = DGsolv_calc(4) - DGsolv_calc(3) - DGsolv_calc(2);

        % answer file
        fprintf(fileID,'%150s %50s %100s %150s %100s %150s %100s %150s %100s %30.0f %30.0f %30.5f %30.5f %30.5f %30.5f %30.5f %30.5f %30.5f %30.5f\n',...
            molecule_list{1,1}.name,molecule_list{1,1}.cas,molecule_list{1,1}.family,...
            molecule_list{2,1}.name,molecule_list{2,1}.family,...
            molecule_list{3,1}.name,molecule_list{3,1}.family,...
            molecule_list{4,1}.name,molecule_list{4,1}.family,...
            Tk,Pbar,DGsolv_exp_1,DGsolv_calc(2),...
            DGsolv_exp_2,DGsolv_calc(3),...
            DGsolv_exp_TS,DGsolv_calc(4),...
            DDG_exp,DDG_calc);
        Esolv_res.DDG_exp(i) = DDG_exp;
        Esolv_res.DDG_calc(i) = DDG_calc;
    else
        Esolv_res.DDG_exp(i) = 8888;
        Esolv_res.DDG_calc(i) = 8888;
    end

end

delete(f);
fclose(fileID);

end