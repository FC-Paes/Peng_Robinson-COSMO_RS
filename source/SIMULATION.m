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

function Esolv_res = SIMULATION(Esolv,data,par)

%--------------------------------------------------------------------------
% file head (answer file)
fileID = fopen('solvation-data\solv-energy-res.dat','w');
fprintf(fileID,'%150s %50s %150s %150s %50s %150s %30s %30s %30s %30s %30s %30s\n',...
    'Name_solvent','CAS_solvent','Family_solvent',...
    'Name_solute','CAS_solute','Family_solute',...
    'BAC','Type_BAC','T[K]','P[bar]',...
    'DGsolv_exp[kcal/mol]','DGsolv_calc[kcal/mol]');

% Number of experimental points to be simulated
Npoints = size(Esolv,1);

% wait bar to follow the progress of the calculation
msg = strcat('calculation progress:',num2str(0),' out of ',num2str(Npoints));
f = waitbar(0,msg,'Name','Calculation of solvation free energies');

%--------------------------------------------------------------------------
% Main loop
Esolv_res = Esolv;

% Initialization of the COSMO-RS interaction matrix at 298 K
[par.Em_298, par.Ehb_298] = INTERACTION_MATRIX_298(par);

for i = 1:Npoints

    %-------------------------
    % update of the wait bar
    msg = strcat('calculation progress:',num2str(i),'-->',num2str(Npoints));
    waitbar(i/Npoints,f,sprintf(msg))

    %-------------------------
    % Point information (solvent, solute, temperature, and pressure)
    solvent = string(Esolv.solvent(i));
    solute = string(Esolv.solute(i));
    Tk = Esolv.T_K_(i);
    isHead = ismember('P_bar_', Esolv.Properties.VariableNames);
    if isHead
        Pbar = Esolv.P_bar_(i);
    else
        Pbar = 0; % in this case we are going to use the saturation pressure
    end
    DGsolv_exp = Esolv.DGsolv_kcal_mol_(i);

    %-------------------------
    % Check if both solvent and solute are present in the select databases
    % (both PR parameters and sigma-profiles must be available or able to be
    % predicted through group contribution)

    %===========
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
    % --> solute
    stat_solute = 0;

    %   - check PR inputs
    if par.EoS_par.solute == "Exp data"
        % check if PR EoS inputs are available in the experimental database
        stat_solute_crit = ismember(solute,data.exp.Name_simulis_);
    else
        % Ohterwise we check the group decomposition
        stat_solute_crit = ismember(solute,data.gc.Name_simulis_);
    end

    %   - check sigma profiles
    if par.EoS_par.solute == "Exp data"
        % check if sigma-profiles are available in the quantum-based database
        stat_solute_sp = 0;
        filename = strcat('sigma-profiles\',solute,'.dat');
        if isfile(filename)
            stat_solute_sp = 1;
        end
    else
        % Ohterwise we check the group decomposition
        stat_solute_sp = ismember(solute,data.gc.Name_simulis_);
    end

    % Can the solute be used?
    % -- stat_solute = 0 (No)
    % -- stat_solute = 1 (yes)
    if stat_solute_sp == 1 && stat_solute_crit == 1
        stat_solute = 1;
    end

    %---------------------------------------------
    % main loop
    if stat_solvent == 1 && stat_solute == 1

        %--------------------------
        IDtype = "name";

        % solvent
        molecule_list{1,1} = DATA_MOL(solvent,IDtype,data,par.EoS_par.solvent,par.psigma_ipt.solvent,par);

        % solute
        molecule_list{2,1} = DATA_MOL(solute,IDtype,data,par.EoS_par.solute,par.psigma_ipt.solute,par);

        [BAC,Type] = GET_BAC(molecule_list{1,1}.AssociationCode,molecule_list{2,1}.AssociationCode);
        %--------------------------

        % -- composition
        % (1) = solvent
        % (2) = solute

        z=[1;0];

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
        % just in case... (if the aswer is a negative pressure)
        if Pcalc <= 0
            Pcalc = 1; % bar
        end

        %--------------------------
        % calculation of solvation free energies
        phase = 1; % liquid phase
        DGsolv_calc = SOLVATION(Tk,Pcalc,z,phase,molecule_list,par);

        % answer file
        fprintf(fileID,'%150s %50s %150s %150s %50s %150s %30.0f %30.0f %30.5f %30.5f %30.5f %30.5f\n',...
            molecule_list{1,1}.name,molecule_list{1,1}.cas,molecule_list{1,1}.family,...
            molecule_list{2,1}.name,molecule_list{2,1}.cas,molecule_list{2,1}.family,...
            BAC,Type,...
            Tk,Pbar,DGsolv_exp,DGsolv_calc(2));
        Esolv_res.DGsolv_calc_kcal_mol_(i) = DGsolv_calc(2);
    end
    
end
delete(f);
fclose(fileID);

end