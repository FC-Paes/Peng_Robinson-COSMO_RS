%--------------------------------------------------------------------------
% @author: Francisco Carlos Paes, January 2025
% 
% Equipe Thermodynamique et Energie (ThermE)
% Laboratoire Réactions et Génie des Procédés (LRGP)
% UMR 7274 CNRS - Université de Lorraine
%
%==========================================================================
%
% This function optimizes the parametrization of the GC methods for free
% radicals
%
%==========================================================================
% INPUTS
% -- Esolv = data containing Esolv calculated from COSMOTherm
% -- data = structure containing pure compound inputs
% -- par = structure containing the parametrization of tc-PR/COSMO-RS
% -- method = COSMO or GCM (group contribution method) [CHAR 1x1]
% -- radical_for_opt = type of free radical or transition states used 
%                      considered in the optimization
%
% OUTPUTS
% -- parNew = structure containing the new optimized parametrization of the
%             groups considered
% -- results_table = statistics of the optimization
%==========================================================================

function [parNEW,results_table] = OPTIMIZATION_GC(Esolv,data,par,radical_for_opt)

    % groups for the optimisation
    % (see the position of groups in the excel input file)
    pos = 138:176; % all groups
    switch radical_for_opt
        case "H-atom"
            pos = 140;
        case "Acetyllenic-radical"
            pos = 141:142;
        case "Alkoxy-radical"
            pos = 143:148;
        case "Benzyl-radical"
            pos = 149:150;
        case "Carbonyl-radical"
            pos = 151:158;
        case "Peroxy-radical"
            pos = 159:165;
        case "Phenyl-radical"
            pos = 166;
        case "Primary-radical"
            pos = 167:173;
        case "Secondary-radical"
            pos = 174:177;
        case "Tertiary-radical"
            pos = 178;
        case "Vinyl-radical"
            pos = 179:185;
        case "Transition-state"
            pos = 186:205;
    end

    %-------------------------------------
    % Initial guess (current parameters)
    c0(:,1) = par.GCM.ac(pos);
    c0(:,2) = par.GCM.b(pos);
    c0(:,3) = par.GCM.a0(pos);

    %-------------------------------------
    % Optimization (fine tuned parameters based on solvation free energies)
    options = optimset('Display','iter','PlotFcns',@optimplotfval);
    fun = @(c)optunivpar(c,Esolv,data,par,pos);
    c = fminsearch(fun,c0,options);

    %-------------------------------------
    % update of the new parametriztion with the fine-tuned group
    % contribution parameters
    parNEW = par;
    parNEW.GCM.ac(pos) = c(:,1);
    parNEW.GCM.b(pos) = c(:,2);
    parNEW.GCM.a0(pos) = c(:,3);

    %-------------------------------------
    % print results on the command window
    all_groups = par.GCM.groups_UNIFAC;
    groups = all_groups(pos);
    parameters_ac = c(:,1);
    parameters_b = c(:,2);
    parameters_a0 = c(:,3);
    results_table = table(groups,parameters_ac,parameters_b,parameters_a0);

%==========================================================================
% OBJECTIVE FUNCTION
%==========================================================================
    function obj = optunivpar(c,Esolv,data,par,pos)

        %-------------------------------------
        % Save intermadiate results in the 'opt_values.txt' file
        [fid,msg] = fopen('opt_values.txt','wt');
        assert(fid>=3,msg);
        fprintf(fid,'%20.5f\n',c);
        fclose(fid);

        %-------------------------------------
        % Group contribution parameters
        parINT = par;
        parINT.GCM.ac(pos) = c(:,1);
        parINT.GCM.b(pos) = c(:,2);
        parINT.GCM.a0(pos) = c(:,3);

        %-------------------------------------
        % Simulation
        Esolv_res = SIMULATION(Esolv,data,parINT);

        %-------------------------------------
        % Simulation
        data_raw(:,1) = Esolv_res.DGsolv_kcal_mol_;
        data_raw(:,2) = Esolv_res.DGsolv_calc_kcal_mol_;

        %-------------------------------------
        % Correction (elimination of NaN results)
        data = data_raw(data_raw(:,2)~=8888,:);
        exp = data(:,1);
        calc = data(:,2);      

        %-------------------------------------
        % Objective function calculation
        obj = mean(abs(exp - calc));

    end
%==========================================================================
end