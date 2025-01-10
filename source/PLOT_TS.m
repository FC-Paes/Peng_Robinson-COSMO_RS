function PLOT_TS(data_raw)

unit = 'kcal/mol';
%==========================================================================
% data
%==========================================================================
% Creating a array with all results array
% % some points in which the calculation was not possible to be done, we got
% % as result "NaN". These data need to be removed from the results before
% % plotting it !
data_all = data_raw(data_raw.DGsolv_calc_kcal_mol_~=8888,:);

% Data divided by type of free radical
Types_radicals = unique(data_all.Family_solute);

data.ROR  = [];
data.ROOR  = [];

for i = 1:size(Types_radicals,1)
    type = string(Types_radicals(i));
    data_type = data_all((data_all.Family_solute == type),:);
    if type == "RORp-Transition-state" || type == "RORs-Transition-state" || type == "RORt-Transition-state"
        data.ROR = [data.ROR;data_type];
    elseif type == "ROORp-Transition-state" || type == "ROORs-Transition-state" || type == "ROORt-Transition-state"
        data.ROOR = [data.ROOR;data_type];
    end
end

AAD.ROR = mean(abs(data.ROR.DGsolv_exp_kcal_mol_ - data.ROR.DGsolv_calc_kcal_mol_));
AAD.ROOR = mean(abs(data.ROOR.DGsolv_exp_kcal_mol_ - data.ROOR.DGsolv_calc_kcal_mol_));

%==========================================================================
% Plots
%==========================================================================

% Define formatting specification
formatSpec = '%.2f';

figure()

%--------------------------------------------------------------------------
% Type = 2 (C-centered radicals - Aryl and benzyl)

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

% all data
trueData = data_all.DGsolv_exp_kcal_mol_;
calculatedData = data_all.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
X = data_all.DGsolv_exp_kcal_mol_;
Y = data_all.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(1), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(1));
hold on

% Add bisector line
eps = 0.3;
min_val = min([min(trueData), min(calculatedData)]);
max_val = max([max(trueData), max(calculatedData)]);
x_limits = [floor(min_val - eps*(max_val - min_val)), ceil(max_val + eps*(max_val - min_val))];
y_limits = x_limits;
plot(x_limits, y_limits, '--', 'LineWidth', 2.0,Color="r");
hold on

% Add labels
xlabel('COSMOTherm database', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('tc-PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf(b) transition-states'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' ',unit,')');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    'Data',...
    'Bisector line',...
    'Location','northwest','EdgeColor','none','fontsize',12)

% Set plot limits and grid
xlim("padded");
ylim("padded");

ax = gca;
ax.FontSize = 14;
ax.XAxis.TickLabelFormat = '%.1f';
ax.YAxis.TickLabelFormat = '%.1f';

set(gca, 'Box', 'on', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');

hold off


%--------------------------------------------------------------------------

figure()

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

% all data
data_type = [data.ROR;data.ROOR];
trueData = data_type.DGsolv_exp_kcal_mol_;
calculatedData = data_type.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
% -- ROR
X = data.ROR.DGsolv_exp_kcal_mol_;
Y = data.ROR.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(2), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(2));
hold on

% -- ROOR
X = data.ROOR.DGsolv_exp_kcal_mol_;
Y = data.ROOR.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(1), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(1));
hold on

% Add bisector line
eps = 0.3;
min_val = min([min(trueData), min(calculatedData)]);
max_val = max([max(trueData), max(calculatedData)]);
x_limits = [floor(min_val - eps*(max_val - min_val)), ceil(max_val + eps*(max_val - min_val))];
y_limits = x_limits;
plot(x_limits, y_limits, '--', 'LineWidth', 2.0,Color="r");
hold on

% Add labels
xlabel('COSMOTherm database', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('tc-PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf H-abstraction transition-states'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' kcal/mol)');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    strcat('Perfomerd by alkoxyl radical (MUE = ',num2str(AAD.ROR,formatSpec),')'),...
    strcat('Perfomerd by peroxyl radical (MUE = ',num2str(AAD.ROOR,formatSpec),')'),...
        'Bisector line',...
    'Location','northwest','EdgeColor','none','fontsize',11)

% Set plot limits and grid
xlim("padded");
ylim("padded");

ax = gca;
ax.FontSize = 14;
ax.XAxis.TickLabelFormat = '%.1f';
ax.YAxis.TickLabelFormat = '%.1f';

set(gca, 'Box', 'on', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');

hold off;
