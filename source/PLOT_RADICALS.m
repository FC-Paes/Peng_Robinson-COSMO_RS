function PLOT_RADICALS(data_raw)

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

data.alkyl  = [];
data.vinyl  = [];
data.acetylenic  = [];
data.benzyl  = [];
data.aryl  = [];
data.hydroxy  = [];
data.alkoxyl  = [];
data.peroxy  = [];
data.acyl  = [];
data.acyloxy  = [];


for i = 1:size(Types_radicals,1)
    type = string(Types_radicals(i));
    data_type = data_all((data_all.Family_solute == type),:);

    % alkyl vinyl, and acetylenic
    if type == "Primary-radical" || type == "Secondary-radical" || type == "Tertiary-radical"
        data.alkyl = [data.alkyl;data_type];
        AAD.alkyl = mean(abs(data.alkyl.DGsolv_exp_kcal_mol_ - data.alkyl.DGsolv_calc_kcal_mol_));
    elseif type == "Vinyl-radical" 
        data.vinyl = [data.vinyl;data_type];
        AAD.vinyl = mean(abs(data.vinyl.DGsolv_exp_kcal_mol_ - data.vinyl.DGsolv_calc_kcal_mol_));
    elseif type == "Acetylenic-radical"
        data.acetylenic = [data.acetylenic;data_type];
        AAD.acetylenic = mean(abs(data.acetylenic.DGsolv_exp_kcal_mol_ - data.acetylenic.DGsolv_calc_kcal_mol_));

    % aryl and benzyl
    elseif type == "Benzyl-radical"
        data.benzyl = [data.benzyl;data_type];
        AAD.benzyl = mean(abs(data.benzyl.DGsolv_exp_kcal_mol_ - data.benzyl.DGsolv_calc_kcal_mol_));
     elseif type == "Aryl-radical"
        data.aryl = [data.aryl;data_type];
        AAD.aryl = mean(abs(data.aryl.DGsolv_exp_kcal_mol_ - data.aryl.DGsolv_calc_kcal_mol_));

    % alkoxy and peroxy
    elseif type == "Hydroxyl-radical"
        data.hydroxy = [data.hydroxy;data_type];
        AAD.hydroxy = mean(abs(data.hydroxy.DGsolv_exp_kcal_mol_ - data.hydroxy.DGsolv_calc_kcal_mol_));
    elseif type == "Alkoxy-radical"
        data.alkoxyl = [data.alkoxyl;data_type];
        AAD.alkoxyl = mean(abs(data.alkoxyl.DGsolv_exp_kcal_mol_ - data.alkoxyl.DGsolv_calc_kcal_mol_));
    elseif type == "Peroxy-radical" 
        data.peroxy = [data.peroxy;data_type];
        AAD.peroxy = mean(abs(data.peroxy.DGsolv_exp_kcal_mol_ - data.peroxy.DGsolv_calc_kcal_mol_));

    % Carbonyl
    elseif type == "Acyl-radical"
        data.acyl = [data.acyl;data_type];
        AAD.acyl = mean(abs(data.acyl.DGsolv_exp_kcal_mol_ - data.acyl.DGsolv_calc_kcal_mol_));   
    elseif type == "Acyloxy-radical"
        data.acyloxy = [data.acyloxy;data_type];
        AAD.acyloxy = mean(abs(data.acyloxy.DGsolv_exp_kcal_mol_ - data.acyloxy.DGsolv_calc_kcal_mol_));   
    end
end

%==========================================================================
% Plots
%==========================================================================

% Define formatting specification
formatSpec = '%.2f';

figure()

%--------------------------------------------------------------------------
% Type = 1 (C-centered radicals - alkyl, vinyl, and acetylenic)

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

subplot(2,2,1)

% all data
data_type = [data.alkyl;data.vinyl;data.acetylenic];
trueData = data_type.DGsolv_exp_kcal_mol_;
calculatedData = data_type.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
% -- alkyl
X = data.alkyl.DGsolv_exp_kcal_mol_;
Y = data.alkyl.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(1), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(1));
hold on

% -- vinyl
X = data.vinyl.DGsolv_exp_kcal_mol_;
Y = data.vinyl.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(2), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(2));

% -- acetylenic
X = data.acetylenic.DGsolv_exp_kcal_mol_;
Y = data.acetylenic.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(3), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(3));

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
ylabel('PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf(a) Carbon-Centered Radicals (aliphatics)'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' ',unit,')');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    strcat('sp³ hybridized (MUE = ',num2str(AAD.alkyl,formatSpec),')'),...
    strcat('sp² hybridized (MUE = ',num2str(AAD.vinyl,formatSpec),')'),...
    strcat('sp  hybridized (MUE = ',num2str(AAD.acetylenic,formatSpec),')'),...
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
% Type = 2 (C-centered radicals - Aryl and benzyl)

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

subplot(2,2,2)

% all data
data_type = [data.aryl;data.benzyl];
trueData = data_type.DGsolv_exp_kcal_mol_;
calculatedData = data_type.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
% -- aryl
X = data.aryl.DGsolv_exp_kcal_mol_;
Y = data.aryl.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(2), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(2));
hold on

% -- benzyl
X = data.benzyl.DGsolv_exp_kcal_mol_;
Y = data.benzyl.DGsolv_calc_kcal_mol_;
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
ylabel('PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf(b) Carbon-Centered Radicals (aromatics)'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' ',unit,')');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    strcat('Aryl radicals (MUE = ',num2str(AAD.aryl,formatSpec),')'),...
    strcat('Benzyl radicals (MUE = ',num2str(AAD.benzyl,formatSpec),')'),...
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
% Type = 3 (Carbonyl)

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

subplot(2,2,3)

% all data
data_type = [data.acyl;data.acyloxy];
trueData = data_type.DGsolv_exp_kcal_mol_;
calculatedData = data_type.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
% -- aciloxy
X = data.acyloxy.DGsolv_exp_kcal_mol_;
Y = data.acyloxy.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(1), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(1));
hold on
% -- acyl
X = data.acyl.DGsolv_exp_kcal_mol_;
Y = data.acyl.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(2), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(2));

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
ylabel('PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf(c) Radicals formed on carbonyl groups'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' ',unit,')');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    strcat("Acyloxyl radicals (MUE = ",num2str(AAD.acyloxy,formatSpec),')'),...
    strcat("Acyl radicals (MUE = ",num2str(AAD.acyl,formatSpec),')'),...
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

% Type = 4 (O-centered radicals - alkoxyl and peroxyl)

marker=['o','s','^','d'];
color=["#0072BD","#D95319","#77AC30","#EDB120"];

subplot(2,2,4)

% all data
data_type = [data.hydroxy;data.alkoxyl;data.peroxy];
trueData = data_type.DGsolv_exp_kcal_mol_;
calculatedData = data_type.DGsolv_calc_kcal_mol_;
AAD_type = mean(abs(trueData - calculatedData));

% Add scatter plots
% -- alkoxy
X = data.alkoxyl.DGsolv_exp_kcal_mol_;
Y = data.alkoxyl.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(1), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(1));
hold on

% -- peroxy
X = data.peroxy.DGsolv_exp_kcal_mol_;
Y = data.peroxy.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(2), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(2));
hold on

% -- hydroxy
X = data.hydroxy.DGsolv_exp_kcal_mol_;
Y = data.hydroxy.DGsolv_calc_kcal_mol_;
scatter(X, Y, 30, marker(3), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, ...
    'MarkerEdgeColor','black', ...
    'MarkerFaceColor',color(3));
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
ylabel('PR/COSMO-RS', 'FontSize', 14, 'FontWeight', 'bold');

% Add title
format1 = ['\fontsize{14}\color{black}\bf(d) Radicals formed on hydroxyl groups'];
format2 = strcat('\fontsize{14}\color{black}\rm','(MUE=',num2str(AAD_type,formatSpec),' ',unit,')');
formattedText = {format1;format2};
title(formattedText)

% Add legend
legend(...
    strcat("alkoxyl radicals (MUE = ",num2str(AAD.alkoxyl,formatSpec),')'),...
    strcat("peroxyl radicals (MUE = ",num2str(AAD.peroxy,formatSpec),')'),...
    strcat("hydroxyl radical (MUE = ",num2str(AAD.hydroxy,formatSpec),')'),...
    "Bisector line",...
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

% Deviations
fprintf('\n=== AAD Values for Free Radicals ===\n');
radicalTypes = fieldnames(AAD);

for i = 1:length(radicalTypes)
    type = radicalTypes{i};
    value = AAD.(type);
    if isnan(value)
        fprintf('%s: No valid data (contains NaN)\n', type);
    else
        fprintf('%s: %.2f kcal/mol\n', type, value);
    end
end
fprintf('=====================================\n');

