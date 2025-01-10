function PLOT_ACTIVATION(data_raw)

%==========================================================================
% data
%==========================================================================
% Creating a array with all results array
% % some points in which the calculation was not possible to be done, we got
% % as result "NaN". These data need to be removed from the results before
% % plotting it !
data_all = data_raw(data_raw.DDGsolv_calc_kcal_mol_~=8888,:);
data_298 = data_all(data_all.T_K_==298,:);

% Data divided by type of free radical
Types_radicals = unique(data_all.Family_TS);

data.ROR  = [];
data.ROOR  = [];

for i = 1:size(Types_radicals,1)
    type = string(Types_radicals(i));
    data_type = data_all((data_all.Family_TS == type),:);
    if type == "RORp-Transition-state" || type == "RORs-Transition-state" || type == "RORt-Transition-state"
        data.ROR = [data.ROR;data_type];
    elseif type == "ROORp-Transition-state" || type == "ROORs-Transition-state" || type == "ROORt-Transition-state"
        data.ROOR = [data.ROOR;data_type];
    end
end

AAD.ROR = mean(abs(data.ROR.DDGsolv_exp_kcal_mol_ - data.ROR.DDGsolv_calc_kcal_mol_));
AAD.ROOR = mean(abs(data.ROOR.DDGsolv_exp_kcal_mol_ - data.ROOR.DDGsolv_calc_kcal_mol_));

%==========================================================================
% Plots
%==========================================================================

% Define formatting specification
formatSpec = '%.2f';

figure()

%--------------------------------------------------------------------------

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

% all data
trueData = data_all.DDGsolv_exp_kcal_mol_;
calculatedData = data_all.DDGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
X = data_all.DDGsolv_exp_kcal_mol_;
Y = data_all.DDGsolv_calc_kcal_mol_;
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
format1 = ['\fontsize{14}\color{black}\bf(a) Solvation Gibbs energies of activation'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' kcal/mol)');
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

hold off;

%--------------------------------------------------------------------------
% Plot error distribution

figure()

% Calculate error (difference between true and calculated values)
errorData = data_all.DDGsolv_calc_kcal_mol_ - data_all.DDGsolv_exp_kcal_mol_;
errorData_298 = data_298.DDGsolv_calc_kcal_mol_ - data_298.DDGsolv_exp_kcal_mol_;

% Calculate standard deviation of the error
std_dev = std(errorData);

% Plot histogram of error distribution
histogram(errorData, 'FaceColor', "#0072BD", 'EdgeColor', 'black', 'FaceAlpha', 0.6);
hold on
histogram(errorData_298, 'FaceColor', "#F28C28", 'EdgeColor', 'black', 'FaceAlpha', 0.6);
hold on

% Add a vertical line
xline(0, '--r', 'LineWidth', 2, 'DisplayName', 'No Error');

% Add labels
xlabel('Error (kcal/mol)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');

% Add title with MAD
format1 = '\fontsize{14}\color{black}\bf(b) Error Distribution';
format2 = strcat('\fontsize{14}\color{black}\rm','(SD=', num2str(std_dev, formatSpec), ' kcal/mol)');
formattedText = {format1; format2};
title(formattedText)

% Add legend
legend('All Temperatures','298 K only','Zero-centered line', 'Location', 'northwest', 'EdgeColor', 'none', 'fontsize', 12);

% Set plot limits and grid
xlim("padded");
grid off

ax = gca;
ax.FontSize = 14;
set(gca, 'Box', 'on', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');

hold off;
