function PLOT_MOLECULES(data_raw)

%==========================================================================
% data
%==========================================================================
% Creating a array with all results array
% % some points in which the calculation was not possible to be done, we got
% % as result "NaN". These data need to be removed from the results before
% % plotting it !

data = data_raw(data_raw.DGsolv_calc_kcal_mol_~=8888,:);

% Average absolute deviation in kcal/mol
AAD = mean(abs(data.DGsolv_exp_kcal_mol_ - data.DGsolv_calc_kcal_mol_));
fprintf('\n=====================================\n');
fprintf('\n Average absolute deviations (AAD) \n');
fprintf('\n');
fprintf('\n=== All data ===\n');
fprintf('AAD = %.4f kcal/mol\n', AAD);

% Data divided by binary association code (BAC)
fprintf('\n');
fprintf('\n=== Binary association code ===\n');
for i = 1:9
    % data
    BAC{i,1}.data = data((data.BAC == i),:);
    % AAD
    BAC{i,1}.AAD = mean(abs(BAC{i,1}.data.DGsolv_exp_kcal_mol_ - BAC{i,1}.data.DGsolv_calc_kcal_mol_));
    % Print the result
    fprintf('BAC %d: AAD = %.4f kcal/mol\n', i, BAC{i,1}.AAD);
end

% calculation of AAD by type of association
fprintf('\n');
fprintf('\n=== Type of association ===\n');
fprintf('* type 1 = no association\n');
fprintf('* type 2 = only self association takes place\n');
fprintf('* type 3 = only cross association takes place\n');
fprintf('* type 4 = self and cross association take place\n');
fprintf(' \n');

% Initialize AAD_type
AAD_type = zeros(4, 1);

% Type 1 calculation
numerator = 0;
denominator = 0;
for i = 1:4
    if ~isnan(BAC{i}.AAD)
        numerator = numerator + size(BAC{i}.data, 1) * BAC{i}.AAD;
        denominator = denominator + size(BAC{i}.data, 1);
    end
end
AAD_type(1, 1) = numerator / denominator;

% Type 2 calculation
if ~isnan(BAC{5}.AAD)
    AAD_type(2, 1) = BAC{5}.AAD;
else
    AAD_type(2, 1) = NaN; % Or another placeholder if all values are NaN
end

% Type 3 calculation
if ~isnan(BAC{6}.AAD)
    AAD_type(3, 1) = BAC{6}.AAD;
else
    AAD_type(3, 1) = NaN; % Or another placeholder
end

% Type 4 calculation
numerator = 0;
denominator = 0;
for i = 7:9
    if ~isnan(BAC{i}.AAD)
        numerator = numerator + size(BAC{i}.data, 1) * BAC{i}.AAD;
        denominator = denominator + size(BAC{i}.data, 1);
    end
end
AAD_type(4, 1) = numerator / denominator;



for i = 1:4
    % Print the result
    fprintf('Type %d: AAD = %.4f kcal/mol\n', i, AAD_type(i,1));
end
fprintf('=====================================\n');


%==========================================================================
% Plots
%==========================================================================

figure()

marker=['o','s','^','d'];
color=["#0096FF","#FFA500","#A9A9A9","#FFEA00"];
legendLabels = ["Type 4: Self and cross-association","Type 3: Only cross-association" ,"Type 2: Only self-association","Type 1: No association"];

for i = 4:-1:1
    data_type = data((data.Type_BAC == i),:);
    trueData = data_type.DGsolv_exp_kcal_mol_;
    calculatedData = data_type.DGsolv_calc_kcal_mol_;
    scatter(trueData, calculatedData, 30, marker(i), 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.6, ...
        'MarkerEdgeColor',color(i), ...
        'MarkerFaceColor',color(i));
    hold on;
end

trueData = data.DGsolv_exp_kcal_mol_;
calculatedData = data.DGsolv_calc_kcal_mol_;

% Plot y=x line for reference
eps = 0.15;
min_val = min([min(trueData), min(calculatedData)]);
max_val = max([max(trueData), max(calculatedData)]);
x_limits = [floor(min_val - eps*(max_val - min_val)), ceil(max_val + eps*(max_val - min_val))];
y_limits = x_limits;
plot(x_limits, y_limits, '--', 'LineWidth', 2.0,Color="r");

% Add labels and title
xlabel('CompSol database', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');
title('Solvation Free Energies [kcal/mol]', 'FontSize', 14, 'FontWeight', 'bold');

% Set plot limits and grid
xlim("padded");
ylim("padded");

% Add the legend
legend(legendLabels, 'FontSize', 10, 'Location', 'northwest', 'Box', 'off');

ax = gca;
ax.FontSize = 14;
ax.XAxis.TickLabelFormat = '%.1f';
ax.YAxis.TickLabelFormat = '%.1f';

set(gca, 'Box', 'on', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');

hold off;

end % function

