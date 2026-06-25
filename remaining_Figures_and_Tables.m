%% Multi-subject data aggregation and adaptation-effect analysis
% this script include
% Figure 2A: Change in natural logarithmic preferred numerosity as a function
% Figure 3A+B: Comparison of distance-dependent numerosity adaptation
% Figure S3A+B: Linear-scale benchmarking of distance-dependent modulation
% Table S2: Voxel-pooled bootstrap test of the condition-difference index (Diff) for distance-dependent adaptation slopes
% Table S3: Robustness of voxel-pooled bootstrap estimates
%  
clear; clc; close all;

%% Data directory
data_dir = 'D:\lrs\all-in-one analysis\params';
mat_files = dir(fullfile(data_dir, '*_merged_data.mat'));
if isempty(mat_files)
    error('No *_merged_data.mat files were found in "%s". Please check the path or file names.', data_dir);
end

subject_list = cell(1, numel(mat_files));
for i = 1:numel(mat_files)
    [~, name, ~] = fileparts(mat_files(i).name);
    underscore_pos = strfind(name, '_merged_data');
    if isempty(underscore_pos)
        error('File name "%s" does not match the pattern "*_merged_data.mat".', mat_files(i).name);
    end
    subject_list{i} = name(1:underscore_pos-1);
end

fprintf('Found %d subject files in "%s":\n', data_dir, numel(subject_list));
disp(subject_list);

%% Define ROI list (fixed row order)
roi_list = {'NTO', 'NPO', 'NPC1', 'NPC2', 'NPC3', 'NF'};

% Fixed order for slope storage (6x2)
roi_order  = roi_list;
cond_order = {'unconnect','connect'};

% Output directory for saving slopes
slope_out_dir = fullfile(data_dir, 'fit_slopes_saved');
if ~exist(slope_out_dir, 'dir'); mkdir(slope_out_dir); end

% Initialize merged data structures
merged_pref = struct();
merged_ve   = struct();

for r = 1:length(roi_list)
    roi = roi_list{r};
    merged_pref.(roi).con     = [];
    merged_pref.(roi).uncon   = [];
    merged_pref.(roi).control = [];
    merged_ve.(roi).control   = [];
end

subject_voxel_counts = struct();

%% Step 1: iterate through all subjects and merge data
for s = 1:length(subject_list)
    subject  = subject_list{s};
    filename = fullfile(data_dir, [subject '_merged_data.mat']);

    if ~exist(filename, 'file')
        fprintf('File does not exist: %s. Skipping this subject.\n', filename);
        continue;
    end

    try
        fprintf('Loading: %s\n', filename);
        loaded_data = load(filename, 'pref', 've');
        pref = loaded_data.pref;
        ve   = loaded_data.ve;

        for r = 1:length(roi_list)
            roi = roi_list{r};

            if isfield(pref, roi) && isfield(ve, roi) && ...
               isfield(pref.(roi), 'control') && isfield(ve.(roi), 'control')

                con_data     = pref.(roi).con;
                uncon_data   = pref.(roi).uncon;
                control_data = pref.(roi).control;
                ve_data      = ve.(roi).control;

                if ~isempty(con_data) && ~isempty(uncon_data) && ~isempty(control_data) && ...
                   numel(con_data) == numel(uncon_data) && ...
                   numel(con_data) == numel(control_data)

                    merged_pref.(roi).con     = [merged_pref.(roi).con,     con_data];
                    merged_pref.(roi).uncon   = [merged_pref.(roi).uncon,   uncon_data];
                    merged_pref.(roi).control = [merged_pref.(roi).control, control_data];
                    merged_ve.(roi).control   = [merged_ve.(roi).control,   ve_data];

                    if ~isfield(subject_voxel_counts, roi)
                        subject_voxel_counts.(roi) = [];
                    end
                    subject_voxel_counts.(roi) = [subject_voxel_counts.(roi), numel(con_data)];

                    fprintf('  %s: %d voxels\n', roi, numel(con_data));
                else
                    warning('Data dimension mismatch or empty data: %s - %s', subject, roi);
                end
            else
                fprintf('Subject %s is missing ROI %s or required fields.\n', subject, roi);
            end
        end

        fprintf('Successfully loaded and merged subject: %s\n', subject);

    catch ME
        warning('Error while loading subject %s: %s', subject, ME.message);
    end
end

% Check whether any data were successfully loaded
total_voxels = 0;
for r = 1:length(roi_list)
    roi = roi_list{r};
    if ~isempty(merged_pref.(roi).control)
        total_voxels = total_voxels + numel(merged_pref.(roi).control);
    end
end

if total_voxels == 0
    error('No data were successfully loaded. Please check the file path and file names.');
end

fprintf('\nData merging completed. Total merged voxels: %d\n', total_voxels);

%% Step 2: adaptation-effect analysis (linear scale, with SEM) [connected & unconnected adaptation]
results_lin = struct();

% Linear-space slope fitting (Figure 6A)
% col1 = unconnect, col2 = connect; rows follow roi_order
slope_lin_6x2     = nan(numel(roi_order), numel(cond_order));
intercept_lin_6x2 = nan(numel(roi_order), numel(cond_order));

for r = 1:length(roi_list)
    roi = roi_list{r};

    if isempty(merged_pref.(roi).control)
        fprintf('[Linear scale] ROI %s has no data. Skipping.\n', roi);
        continue
    end

    control_data = merged_pref.(roi).control;
    con_data     = merged_pref.(roi).con;
    uncon_data   = merged_pref.(roi).uncon;
    ve_data      = merged_ve.(roi).control;

    fprintf('\n[Linear scale] Analyzing ROI: %s (total voxels: %d)\n', roi, numel(control_data));

    % Data quality filtering
    valid_idx = (ve_data > 0.3) & ...
                (control_data >= 1.05) & (control_data <= 7.05) & ...
                (con_data     >= 1.05) & (con_data     <= 7.05) & ...
                (uncon_data   >= 1.05) & (uncon_data   <= 7.05);

    control_data_filtered = control_data(valid_idx);
    con_data_filtered     = con_data(valid_idx);
    uncon_data_filtered   = uncon_data(valid_idx);

    fprintf('[Linear scale] Before filtering: %d voxels, after filtering: %d voxels (retained %.1f%%)\n', ...
            numel(control_data), numel(control_data_filtered), ...
            numel(control_data_filtered)/numel(control_data)*100);

    if isempty(control_data_filtered)
        warning('[Linear scale] No valid data after filtering for ROI %s. Skipping.', roi);
        continue
    end

    control_data = control_data_filtered;
    con_data     = con_data_filtered;
    uncon_data   = uncon_data_filtered;

    % Binning according to control PN ~ 1-7
    bin_indices = zeros(size(control_data));
    bin_indices(control_data >= 1.05   & control_data <= 1.5) = 1;
    bin_indices(control_data > 1.5 & control_data <= 2.5) = 2;
    bin_indices(control_data > 2.5 & control_data <= 3.5) = 3;
    bin_indices(control_data > 3.5 & control_data <= 4.5) = 4;
    bin_indices(control_data > 4.5 & control_data <= 5.5) = 5;
    bin_indices(control_data > 5.5 & control_data <= 6.5) = 6;
    bin_indices(control_data > 6.5 & control_data <= 7.5) = 7;

    con_diff_means   = nan(1, 7);
    uncon_diff_means = nan(1, 7);
    con_diff_sem     = nan(1, 7);
    uncon_diff_sem   = nan(1, 7);
    bin_counts       = zeros(1, 7);

    for bin = 1:7
        idx = (bin_indices == bin);
        bin_counts(bin) = sum(idx);

        if bin_counts(bin) < 2
            if bin_counts(bin) == 1
                warning('[Linear scale] ROI %s, bin %d contains only 1 voxel; SEM cannot be computed.', roi, bin);
            end
            continue
        end

        con_diff   = con_data(idx)   - control_data(idx);
        uncon_diff = uncon_data(idx) - control_data(idx);

        con_diff_means(bin)   = mean(con_diff);
        uncon_diff_means(bin) = mean(uncon_diff);
        con_diff_sem(bin)     = std(con_diff)   / sqrt(numel(con_diff));
        uncon_diff_sem(bin)   = std(uncon_diff) / sqrt(numel(uncon_diff));

        fprintf('[Linear scale]   Bin %d: %d voxels, con change: %.3f±%.3f, uncon change: %.3f±%.3f\n', ...
                bin, bin_counts(bin), ...
                con_diff_means(bin),   con_diff_sem(bin), ...
                uncon_diff_means(bin), uncon_diff_sem(bin));
    end

    % Save linear-scale results
    results_lin.(roi).con_diff_mean   = con_diff_means;
    results_lin.(roi).uncon_diff_mean = uncon_diff_means;
    results_lin.(roi).con_diff_sem    = con_diff_sem;
    results_lin.(roi).uncon_diff_sem  = uncon_diff_sem;
    results_lin.(roi).bin_counts      = bin_counts;
    results_lin.(roi).total_voxels    = numel(control_data);

    % Visualization (linear scale, with SEM)
    figure('Name', ['[Linear] ROI: ' roi], 'Position', [100, 100, 800, 600]);
    hold on

    x = 1:7;

    valid_con_bins   = ~isnan(con_diff_means);
    valid_uncon_bins = ~isnan(uncon_diff_means);

    errorbar(x(valid_con_bins), con_diff_means(valid_con_bins), con_diff_sem(valid_con_bins), ...
             'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
             'LineWidth', 1.0, 'CapSize', 8, 'DisplayName', 'Connected Adaptation');

    errorbar(x(valid_uncon_bins), uncon_diff_means(valid_uncon_bins), uncon_diff_sem(valid_uncon_bins), ...
             'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
             'LineWidth', 1.0, 'CapSize', 8, 'DisplayName', 'Unconnected Adaptation');

    % Linear fit - con (save slope/intercept to slope_lin_6x2)
    if sum(valid_con_bins) > 1
        p_con    = polyfit(x(valid_con_bins), con_diff_means(valid_con_bins), 1);

        slope_lin_6x2(r, 2)     = p_con(1);
        intercept_lin_6x2(r, 2) = p_con(2);

        yfit_con = polyval(p_con, x(valid_con_bins));

        ymean_con  = mean(con_diff_means(valid_con_bins));
        ss_tot_con = sum((con_diff_means(valid_con_bins) - ymean_con).^2);
        ss_res_con = sum((con_diff_means(valid_con_bins) - yfit_con).^2);
        r2_con     = 1 - (ss_res_con / ss_tot_con);

        b_con = p_con(2);
        if b_con >= 0
            fit_label_con = sprintf('Con fit: y=%.3fx+%.3f (R^2=%.3f)', p_con(1), b_con, r2_con);
        else
            fit_label_con = sprintf('Con fit: y=%.3fx-%.3f (R^2=%.3f)', p_con(1), abs(b_con), r2_con);
        end

        plot(x(valid_con_bins), yfit_con, 'r-', 'LineWidth', 2, ...
             'DisplayName', fit_label_con);

        fprintf('[Linear scale]   Con linear fit R^2 = %.3f\n', r2_con);
    end

    % Linear fit - uncon (save slope/intercept to slope_lin_6x2)
    if sum(valid_uncon_bins) > 1
        p_uncon  = polyfit(x(valid_uncon_bins), uncon_diff_means(valid_uncon_bins), 1);

        slope_lin_6x2(r, 1)     = p_uncon(1);
        intercept_lin_6x2(r, 1) = p_uncon(2);

        yfit_uncon = polyval(p_uncon, x(valid_uncon_bins));

        ymean_uncon  = mean(uncon_diff_means(valid_uncon_bins));
        ss_tot_uncon = sum((uncon_diff_means(valid_uncon_bins) - ymean_uncon).^2);
        ss_res_uncon = sum((uncon_diff_means(valid_uncon_bins) - yfit_uncon).^2);
        r2_uncon     = 1 - (ss_res_uncon / ss_tot_uncon);

        b_uncon = p_uncon(2);
        if b_uncon >= 0
            fit_label_uncon = sprintf('Uncon fit: y=%.3fx+%.3f (R^2=%.3f)', p_uncon(1), b_uncon, r2_uncon);
        else
            fit_label_uncon = sprintf('Uncon fit: y=%.3fx-%.3f (R^2=%.3f)', p_uncon(1), abs(b_uncon), r2_uncon);
        end

        plot(x(valid_uncon_bins), yfit_uncon, 'b-', 'LineWidth', 2, ...
             'DisplayName', fit_label_uncon);

        fprintf('[Linear scale]   Uncon linear fit R^2 = %.3f\n', r2_uncon);
    end

    title(sprintf('[Linear] ROI: %s (n=%d)', roi, numel(control_data)), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Preferred Numerosity (Control Condition)', 'FontSize', 12);
    ylabel('Change in preferred numerosity', 'FontSize', 12);
    legend('Location', 'best');
    grid on; box on
    ylim([-5, 5]);
    set(gca, 'XTick', 1:7, 'XTickLabel', {'1','2','3','4','5','6','7'});
    yline(0, 'k--', 'LineWidth', 1);

    % Annotate the number of voxels in each bin
    for bin = 1:7
        if ~isnan(con_diff_means(bin)) && bin_counts(bin) > 0
            text(bin, con_diff_means(bin)+0.2, sprintf('n=%d', bin_counts(bin)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
        end
    end

    hold off
    fprintf('[Linear scale] Figure created. Save manually if needed.\n');
end

%% Optional: save linear-space slopes (6x2)
% out_lin_mat = fullfile(slope_out_dir, 'slopes_linSpace_6x2.mat');
% out_lin_csv = fullfile(slope_out_dir, 'slopes_linSpace_6x2.csv');
% 
% save(out_lin_mat, 'slope_lin_6x2', 'intercept_lin_6x2', 'roi_order', 'cond_order');
% Tlin = array2table(slope_lin_6x2, 'VariableNames', cond_order, 'RowNames', roi_order);
% writetable(Tlin, out_lin_csv, 'WriteRowNames', true);
% 
% fprintf('\n[Saved] Linear-space slopes:\n  MAT: %s\n  CSV: %s\n', out_lin_mat, out_lin_csv);
% disp(Tlin);

%% Step 3: summary of results (linear scale, optional)
% fprintf('\n=== [Linear scale] Analysis completed ===\n');
% fprintf('ROI summary (linear scale):\n');
% for r = 1:length(roi_list)
%     roi = roi_list{r};
%     if isfield(results_lin, roi)
%         fprintf('  %s: %d valid voxels\n', roi, results_lin.(roi).total_voxels);
%     else
%         fprintf('  %s: no valid data\n', roi);
%     end
% end
% 
% fprintf('\n=== Detailed summary (linear scale, Mean±SEM) ===\n');
% for r = 1:length(roi_list)
%     roi = roi_list{r};
%     if isfield(results_lin, roi)
%         fprintf('\n%s:\n', roi);
%         fprintf('Bin\tCon_Mean±SEM\t\tUncon_Mean±SEM\t\tCount\n');
%         for bin = 1:7
%             if ~isnan(results_lin.(roi).con_diff_mean(bin))
%                 fprintf('%d\t%.3f±%.3f\t\t%.3f±%.3f\t\t%d\n', bin, ...
%                         results_lin.(roi).con_diff_mean(bin), results_lin.(roi).con_diff_sem(bin), ...
%                         results_lin.(roi).uncon_diff_mean(bin), results_lin.(roi).uncon_diff_sem(bin), ...
%                         results_lin.(roi).bin_counts(bin));
%             else
%                 fprintf('%d\tN/A\t\t\tN/A\t\t\t%d\n', bin, results_lin.(roi).bin_counts(bin));
%             end
%         end
%     end
% end
% 
% fprintf('\n[Linear scale] All analyses completed.\n');

%% linear-space three-condition fit visualization (Uncon / Con / High numerosity) -- Figure 6B
% - CSV loading matches Step 5: one CSV per ROI; use the first 7 rows, column 2 as y; x = 1..7
% - Coordinate space matches Step 2: x = 1..7 (control PN bins)
% - Plot style is aligned as closely as possible with the log-scale plotting style in Step 3
% - Revised: fitted lines are always drawn across x = 1..7, extending to bin 7
if ~exist('csv_dir','var')
    csv_dir = 'D:\lrs\all-in-one analysis\high adaptation linear';
end

% Color settings
uncon_color = [0.07, 0.62, 1.00];
con_color   = [1.00, 0.00, 0.00];
gray_high   = [0.50, 0.50, 0.50];

x_lin = 1:7;

fprintf('\n===== Newly added: linear-space three-condition fit plots (Uncon / Con / High) =====\n');

for r = 1:length(roi_list)
    roi = roi_list{r};

    % Check required fields
    if ~isfield(results_lin, roi) || ~isfield(results_lin.(roi), 'con_diff_mean') || isempty(results_lin.(roi).con_diff_mean)
        fprintf('[New - linear three-condition] ROI %s has no results_lin data. Skipping.\n', roi);
        continue;
    end

    % Uncon / Con from Step 2
    yC   = results_lin.(roi).con_diff_mean(:)';
    yU   = results_lin.(roi).uncon_diff_mean(:)';
    semC = results_lin.(roi).con_diff_sem(:)';
    semU = results_lin.(roi).uncon_diff_sem(:)';

    validC = isfinite(yC) & isfinite(semC);
    validU = isfinite(yU) & isfinite(semU);

    % High numerosity from CSV
    yH = nan(1,7);
    csv_file = fullfile(csv_dir, [roi '.csv']);
    if isfile(csv_file)
        try
            Tcsv = readtable(csv_file);
            Mcsv = table2array(Tcsv);
        catch
            Mcsv = readmatrix(csv_file);
        end

        if size(Mcsv,2) >= 2
            Mcsv = Mcsv(:,1:2);
            Mcsv = Mcsv(~any(isnan(Mcsv),2),:);
            if size(Mcsv,1) >= 7
                yH = Mcsv(1:7,2)';
            else
                warning('[New - linear three-condition] %s has fewer than 7 valid rows. High condition skipped for this ROI.', csv_file);
            end
        else
            warning('[New - linear three-condition] %s has fewer than 2 columns. High condition skipped for this ROI.', csv_file);
        end
    else
        warning('[New - linear three-condition] Cannot find %s. High condition skipped for this ROI.', csv_file);
    end
    validH = isfinite(yH);

    % Linear fits; fitted lines always drawn across x = 1..7
    pU = [NaN NaN]; yhatU_all = nan(size(x_lin));
    if sum(validU) > 1
        xU_fit = x_lin(validU);
        yU_fit = yU(validU);
        pU = polyfit(xU_fit, yU_fit, 1);
        yhatU_all = polyval(pU, x_lin);
    end

    pC = [NaN NaN]; yhatC_all = nan(size(x_lin));
    if sum(validC) > 1
        xC_fit = x_lin(validC);
        yC_fit = yC(validC);
        pC = polyfit(xC_fit, yC_fit, 1);
        yhatC_all = polyval(pC, x_lin);
    end

    pH = [NaN NaN]; yhatH_all = nan(size(x_lin));
    if sum(validH) > 1
        xH_fit = x_lin(validH);
        yH_fit = yH(validH);
        pH = polyfit(xH_fit, yH_fit, 1);
        yhatH_all = polyval(pH, x_lin);
    end

    % Plot
    figure('Name', ['[Linear three-condition] ROI: ' roi], 'Position', [120, 120, 820, 560]);
    hold on;

    if any(validU)
        errorbar(x_lin(validU), yU(validU), semU(validU), ...
            'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', uncon_color, ...
            'MarkerEdgeColor', uncon_color, ...
            'Color', uncon_color, ...
            'LineWidth', 1.5, 'CapSize', 8, ...
            'HandleVisibility', 'off');
    end

    if any(validC)
        errorbar(x_lin(validC), yC(validC), semC(validC), ...
            'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', con_color, ...
            'MarkerEdgeColor', con_color, ...
            'Color', con_color, ...
            'LineWidth', 1.5, 'CapSize', 8, ...
            'HandleVisibility', 'off');
    end

    if any(validH)
        plot(x_lin(validH), yH(validH), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', gray_high, ...
            'MarkerEdgeColor', gray_high, ...
            'Color', gray_high, ...
            'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end

    if all(isfinite(pU))
        plot(x_lin, yhatU_all, '-', 'Color', uncon_color, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    if all(isfinite(pC))
        plot(x_lin, yhatC_all, '-', 'Color', con_color,   'LineWidth', 2, 'HandleVisibility', 'off');
    end
    if all(isfinite(pH))
        plot(x_lin, yhatH_all, '-', 'Color', gray_high,   'LineWidth', 2, 'HandleVisibility', 'off');
    end

    yline(0, 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

    grid off;
    ax = gca;
    ax.Box       = 'off';
    ax.TickDir   = 'in';
    ax.LineWidth = 3.0;
    ax.FontSize  = 16;

    xlim([0.8, 7.2]);
    xticks(1:7);
    xticklabels({'1','2','3','4','5','6','7'});

    ylim([-5, 5]);
    yticks(-5:2:5);

    xl = ax.XLabel; xl.Units = 'normalized'; xl.Position(2) = xl.Position(2) - 0.03;
    yl = ax.YLabel; yl.Units = 'normalized'; yl.Position(1) = yl.Position(1) - 0.03;

    h_leg_uncon = plot(nan, nan, 'o', 'MarkerSize', 14, ...
        'MarkerFaceColor', uncon_color, 'MarkerEdgeColor', uncon_color, ...
        'LineStyle', 'none', 'DisplayName', 'Unconnected Adaptation Condition');

    h_leg_con = plot(nan, nan, 'o', 'MarkerSize', 14, ...
        'MarkerFaceColor', con_color, 'MarkerEdgeColor', con_color, ...
        'LineStyle', 'none', 'DisplayName', 'Connected Adaptation Condition');

    h_leg_high = plot(nan, nan, 'o', 'MarkerSize', 14, ...
        'MarkerFaceColor', gray_high, 'MarkerEdgeColor', gray_high, ...
        'LineStyle', 'none', 'DisplayName', 'High numerosity Condition');

    legend([h_leg_uncon, h_leg_con, h_leg_high], 'Location', 'best');

    hold off;

    fprintf('[New - linear three-condition] ROI=%s done. (High CSV: %s)\n', roi, csv_file);
end

fprintf('===== Newly added: linear-space three-condition fit plots completed =====\n');

%% Step 4: adaptation-effect analysis (log scale: x = logPN - log40) -- Figure 4A
adapter_numerosity = 40;
pn_centers         = 1:7;
x_positions_log    = log(pn_centers) - log(adapter_numerosity);
difference_centers = adapter_numerosity - pn_centers;

% Unified x-axis limits
x_min    = min(x_positions_log);
x_max    = max(x_positions_log);
x_range  = x_max - x_min;
x_margin = 0.05 * x_range;
x_limits = [x_min - x_margin, x_max + x_margin];

% x-axis labels displayed as 40-PN
x_labels = arrayfun(@(d) sprintf('%d', d), difference_centers, 'UniformOutput', false);

results = struct();

% Condition colors: row 1 = Unconnected (blue), row 2 = Connected (red)
cond_colors = [0.07, 0.62, 1.00;
               1.00, 0.00, 0.00];
uncon_color = cond_colors(1,:);
con_color   = cond_colors(2,:);
gray        = [0.50 0.50 0.50];

% Log-space slope storage container (6x2)
slope_log_6x2     = nan(numel(roi_order), numel(cond_order));
intercept_log_6x2 = nan(numel(roi_order), numel(cond_order));

for r = 1:length(roi_list)
    roi = roi_list{r};

    if ~isfield(merged_pref, roi) || isempty(merged_pref.(roi).control)
        fprintf('[Log scale] ROI %s has no data. Skipping.\n', roi);
        continue;
    end

    control_data = merged_pref.(roi).control;
    con_data     = merged_pref.(roi).con;
    uncon_data   = merged_pref.(roi).uncon;
    ve_data      = merged_ve.(roi).control;

    fprintf('\n[Log scale] Analyzing ROI: %s (total voxels: %d)\n', roi, numel(control_data));

    % Data quality filtering
    valid_idx = (ve_data > 0.3) & ...
                (control_data >= 1.05) & (control_data <= 7.05) & ...
                (con_data     >= 1.05) & (con_data     <= 7.05) & ...
                (uncon_data   >= 1.05) & (uncon_data   <= 7.05);

    control_data_filtered = control_data(valid_idx);
    con_data_filtered     = con_data(valid_idx);
    uncon_data_filtered   = uncon_data(valid_idx);

    total_voxels_before = numel(control_data);
    total_voxels_after  = numel(control_data_filtered);

    fprintf('[Log scale] Before filtering: %d voxels, after filtering: %d voxels (retained %.1f%%)\n', ...
            total_voxels_before, total_voxels_after, ...
            total_voxels_after / total_voxels_before * 100);

    if isempty(control_data_filtered)
        warning('[Log scale] No valid data after filtering for ROI %s. Skipping.', roi);
        continue;
    end

    control_data = control_data_filtered;
    con_data     = con_data_filtered;
    uncon_data   = uncon_data_filtered;

    % Binning according to control PN ~ 1-7
    bin_indices = zeros(size(control_data));
    bin_indices(control_data >= 1.05  & control_data <= 1.5) = 1;
    bin_indices(control_data > 1.5 & control_data <= 2.5)   = 2;
    bin_indices(control_data > 2.5 & control_data <= 3.5)   = 3;
    bin_indices(control_data > 3.5 & control_data <= 4.5)   = 4;
    bin_indices(control_data > 4.5 & control_data <= 5.5)   = 5;
    bin_indices(control_data > 5.5 & control_data <= 6.5)   = 6;
    bin_indices(control_data > 6.5 & control_data <= 7.5)   = 7;

    con_diff_log_means   = nan(1, 7);
    uncon_diff_log_means = nan(1, 7);
    con_diff_log_sem     = nan(1, 7);
    uncon_diff_log_sem   = nan(1, 7);
    bin_counts           = zeros(1, 7);

    for bin = 1:7
        idx = (bin_indices == bin);
        bin_counts(bin) = sum(idx);

        if bin_counts(bin) < 2
            if bin_counts(bin) == 1
                warning('[Log scale] ROI %s, bin %d contains only 1 voxel; SEM cannot be computed.', roi, bin);
            end
            continue;
        end

        con_diff_log   = log(con_data(idx))   - log(control_data(idx));
        uncon_diff_log = log(uncon_data(idx)) - log(control_data(idx));

        con_diff_log_means(bin)   = mean(con_diff_log);
        uncon_diff_log_means(bin) = mean(uncon_diff_log);
        con_diff_log_sem(bin)     = std(con_diff_log)   / sqrt(numel(con_diff_log));
        uncon_diff_log_sem(bin)   = std(uncon_diff_log) / sqrt(numel(uncon_diff_log));

        fprintf('[Log scale]   Bin %d (PN≈%d, 40-PN=%d): %d voxels, con log-change: %.3f±%.3f, uncon log-change: %.3f±%.3f\n', ...
                bin, pn_centers(bin), difference_centers(bin), bin_counts(bin), ...
                con_diff_log_means(bin),   con_diff_log_sem(bin), ...
                uncon_diff_log_means(bin), uncon_diff_log_sem(bin));
    end

    results.(roi).con_diff_log_mean   = con_diff_log_means;
    results.(roi).uncon_diff_log_mean = uncon_diff_log_means;
    results.(roi).con_diff_log_sem    = con_diff_log_sem;
    results.(roi).uncon_diff_log_sem  = uncon_diff_log_sem;
    results.(roi).bin_counts          = bin_counts;
    results.(roi).total_voxels        = total_voxels_after;

    % Visualization (log scale, x = logPN - log40, with SEM)
    figure('Name', ['[Log] ROI: ' roi], 'Position', [100, 100, 800, 600]);
    hold on;

    valid_con_bins   = ~isnan(con_diff_log_means);
    valid_uncon_bins = ~isnan(uncon_diff_log_means);

    errorbar(x_positions_log(valid_con_bins), ...
             con_diff_log_means(valid_con_bins), ...
             con_diff_log_sem(valid_con_bins), ...
             'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', con_color, ...
             'MarkerEdgeColor', con_color, ...
             'Color', con_color, ...
             'LineWidth', 1.5, 'CapSize', 8, ...
             'HandleVisibility', 'off');

    errorbar(x_positions_log(valid_uncon_bins), ...
             uncon_diff_log_means(valid_uncon_bins), ...
             uncon_diff_log_sem(valid_uncon_bins), ...
             'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', uncon_color, ...
             'MarkerEdgeColor', uncon_color, ...
             'Color', uncon_color, ...
             'LineWidth', 1.5, 'CapSize', 8, ...
             'HandleVisibility', 'off');

    % Linear fit - con (col2 = connect)
    if sum(valid_con_bins) > 1
        x_fit_con = x_positions_log(valid_con_bins);
        y_fit_con = con_diff_log_means(valid_con_bins);

        p_con    = polyfit(x_fit_con, y_fit_con, 1);
        yhat_con = polyval(p_con, x_fit_con);

        plot(x_fit_con, yhat_con, '-', 'Color', con_color, 'LineWidth', 2, 'HandleVisibility', 'off');

        slope_log_6x2(r, 2)     = p_con(1);
        intercept_log_6x2(r, 2) = p_con(2);

        ymean_con  = mean(y_fit_con);
        ss_tot_con = sum((y_fit_con - ymean_con).^2);
        ss_res_con = sum((y_fit_con - yhat_con).^2);
        r2_con     = 1 - (ss_res_con / ss_tot_con);
        fprintf('[Log scale]   Con linear fit (x = logPN - log40) R^2 = %.3f\n', r2_con);
    end

    % Linear fit - uncon (col1 = unconnect)
    if sum(valid_uncon_bins) > 1
        x_fit_uncon = x_positions_log(valid_uncon_bins);
        y_fit_uncon = uncon_diff_log_means(valid_uncon_bins);

        p_uncon    = polyfit(x_fit_uncon, y_fit_uncon, 1);
        yhat_uncon = polyval(p_uncon, x_fit_uncon);

        plot(x_fit_uncon, yhat_uncon, '-', 'Color', uncon_color, 'LineWidth', 2, 'HandleVisibility', 'off');

        slope_log_6x2(r, 1)     = p_uncon(1);
        intercept_log_6x2(r, 1) = p_uncon(2);

        ymean_uncon  = mean(y_fit_uncon);
        ss_tot_uncon = sum((y_fit_uncon - ymean_uncon).^2);
        ss_res_uncon = sum((y_fit_uncon - yhat_uncon).^2);
        r2_uncon     = 1 - (ss_res_uncon / ss_tot_uncon);
        fprintf('[Log scale]   Uncon linear fit (x = logPN - log40) R^2 = %.3f\n', r2_uncon);
    end

    yline(0, 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

    grid off;

    ax = gca;
    ax.Box       = 'off';
    ax.TickDir   = 'in';
    ax.LineWidth = 3.0;
    ax.FontSize  = 16;

    ylim([-1.2, 1.2]);
    yticks(-1.2:0.4:1.2);

    xlim(x_limits);
    xticks(x_positions_log);
    xticklabels(x_labels);

    xl = ax.XLabel; xl.Units = 'normalized';
    xl.Position(2) = xl.Position(2) - 0.03;

    yl = ax.YLabel; yl.Units = 'normalized';
    yl.Position(1) = yl.Position(1) - 0.03;

    h_leg_uncon = plot(nan, nan, 'o', ...
        'MarkerSize', 14, ...
        'MarkerFaceColor', uncon_color, ...
        'MarkerEdgeColor', uncon_color, ...
        'LineStyle', 'none', ...
        'DisplayName', 'Unconnected Adaptation Condition');

    h_leg_con = plot(nan, nan, 'o', ...
        'MarkerSize', 14, ...
        'MarkerFaceColor', con_color, ...
        'MarkerEdgeColor', con_color, ...
        'LineStyle', 'none', ...
        'DisplayName', 'Connected Adaptation Condition');

    legend([h_leg_uncon, h_leg_con], 'Location', 'best');

    hold off;
    fprintf('[Log scale] Figure created. Save manually if needed.\n');
end

% Save log-space slopes: 6x2 (ROI x condition)
% out_log_mat = fullfile(slope_out_dir, 'slopes_logSpace_6x2.mat');
% out_log_csv = fullfile(slope_out_dir, 'slopes_logSpace_6x2.csv');
% 
% save(out_log_mat, 'slope_log_6x2', 'intercept_log_6x2', 'roi_order', 'cond_order', ...
%     'adapter_numerosity', 'pn_centers', 'x_positions_log');
% 
% Tlog = array2table(slope_log_6x2, 'VariableNames', cond_order, 'RowNames', roi_order);
% writetable(Tlog, out_log_csv, 'WriteRowNames', true);
% 
% fprintf('\n[Saved] Log-space slopes:\n  MAT: %s\n  CSV: %s\n', out_log_mat, out_log_csv);
% disp(Tlog);

%% Bootstrap
% Based on the existing workflow:
% voxel merging -> control-based binning -> bin means -> linear fit
% Perform paired voxel-wise bootstrap resampling with replacement to obtain
% slope distributions and test the con-uncon difference.

%% Table S2
% Adjustable parameters
boot_B        = 5000;
min_bin_vox   = 2;
rng_seed      = 20260123;
make_boot_fig = false;

rng(rng_seed);

% Summary table (one row per ROI)
boot_summary = table('Size',[0 12], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'ROI','N_vox','N_eff','slope_uncon_mean','slope_con_mean','diff_mean','diff_median','CI95_low','CI95_high','p_two','q_BH','sig_q05'});

% Storage containers
boot_slope_uncon_all = cell(numel(roi_list),1);
boot_slope_con_all   = cell(numel(roi_list),1);
boot_diff_all        = cell(numel(roi_list),1);

fprintf('\n================ Bootstrap slope-difference test starts ================\n');
fprintf('B=%d, min_bin_vox=%d, seed=%d\n', boot_B, min_bin_vox, rng_seed);

for r = 1:numel(roi_list)
    roi = roi_list{r};

    if ~isfield(merged_pref, roi) || isempty(merged_pref.(roi).control)
        fprintf('[Bootstrap] ROI=%s has no data. Skipping.\n', roi);
        continue;
    end

    control_data = merged_pref.(roi).control;
    con_data     = merged_pref.(roi).con;
    uncon_data   = merged_pref.(roi).uncon;
    ve_data      = merged_ve.(roi).control;

    valid_idx = (ve_data > 0.3) & ...
                (control_data >= 1.05) & (control_data <= 7.05) & ...
                (con_data     >= 1.05) & (con_data     <= 7.05) & ...
                (uncon_data   >= 1.05) & (uncon_data   <= 7.05);

    control_data = control_data(valid_idx);
    con_data     = con_data(valid_idx);
    uncon_data   = uncon_data(valid_idx);

    N = numel(control_data);
    if N < 10
        fprintf('[Bootstrap] ROI=%s has too few valid voxels (N=%d). Skipping.\n', roi, N);
        continue;
    end

    % Binning according to control PN ~ 1-7
    bin_indices = zeros(size(control_data));
    bin_indices(control_data >= 1.05  & control_data <= 1.5) = 1;
    bin_indices(control_data > 1.5 & control_data <= 2.5)   = 2;
    bin_indices(control_data > 2.5 & control_data <= 3.5)   = 3;
    bin_indices(control_data > 3.5 & control_data <= 4.5)   = 4;
    bin_indices(control_data > 4.5 & control_data <= 5.5)   = 5;
    bin_indices(control_data > 5.5 & control_data <= 6.5)   = 6;
    bin_indices(control_data > 6.5 & control_data <= 7.5)   = 7;

    boot_slope_con   = nan(boot_B,1);
    boot_slope_uncon = nan(boot_B,1);

    % Bootstrap main loop: paired voxel resampling
    for b = 1:boot_B
        boot_idx = randi(N, [N,1]);

        control_bs = control_data(boot_idx);
        con_bs     = con_data(boot_idx);
        uncon_bs   = uncon_data(boot_idx);
        bin_bs     = bin_indices(boot_idx);

        con_m   = nan(1,7);
        uncon_m = nan(1,7);

        for bin = 1:7
            idx = (bin_bs == bin);
            if sum(idx) < min_bin_vox
                continue;
            end
            con_m(bin)   = mean(log(con_bs(idx))   - log(control_bs(idx)));
            uncon_m(bin) = mean(log(uncon_bs(idx)) - log(control_bs(idx)));
        end

        vcon   = ~isnan(con_m);
        vuncon = ~isnan(uncon_m);

        if sum(vcon) > 1
            p = polyfit(x_positions_log(vcon), con_m(vcon), 1);
            boot_slope_con(b) = p(1);
        end
        if sum(vuncon) > 1
            p = polyfit(x_positions_log(vuncon), uncon_m(vuncon), 1);
            boot_slope_uncon(b) = p(1);
        end
    end

    boot_diff = boot_slope_con - boot_slope_uncon;

    % Only use iterations where both slopes are valid
    valid_iter = ~isnan(boot_diff);
    N_eff      = sum(valid_iter);

    if N_eff < 100
        fprintf('[Bootstrap] ROI=%s has too few effective iterations (N_eff=%d); results may be unstable.\n', roi, N_eff);
    end
    if N_eff == 0
        fprintf('[Bootstrap] ROI=%s has no valid iterations. Skipping.\n', roi);
        continue;
    end

    diff_valid        = boot_diff(valid_iter);
    slope_con_valid   = boot_slope_con(valid_iter);
    slope_uncon_valid = boot_slope_uncon(valid_iter);

    CI95 = prctile(diff_valid, [2.5 97.5]);

    % Two-sided p
    p_two = 2 * min(mean(diff_valid >= 0), mean(diff_valid <= 0));

    % Store in results structure
    results.(roi).bootstrap.B                    = boot_B;
    results.(roi).bootstrap.min_bin_vox          = min_bin_vox;
    results.(roi).bootstrap.rng_seed             = rng_seed;
    results.(roi).bootstrap.N_vox                = N;
    results.(roi).bootstrap.N_eff                = N_eff;
    results.(roi).bootstrap.slope_con            = boot_slope_con;
    results.(roi).bootstrap.slope_uncon          = boot_slope_uncon;
    results.(roi).bootstrap.diff_con_minus_uncon = boot_diff;
    results.(roi).bootstrap.CI95                 = CI95;
    results.(roi).bootstrap.p_two                = p_two;

    boot_summary = [boot_summary; {string(roi), N, N_eff, ...
        mean(slope_uncon_valid), mean(slope_con_valid), ...
        mean(diff_valid), median(diff_valid), ...
        CI95(1), CI95(2), p_two, NaN, NaN}]; %#ok<AGROW>

    boot_slope_uncon_all{r} = boot_slope_uncon;
    boot_slope_con_all{r}   = boot_slope_con;
    boot_diff_all{r}        = boot_diff;

    fprintf('[Bootstrap] ROI=%s | N=%d | N_eff=%d | mean(diff)=%.4f | CI95=[%.4f, %.4f] | p=%.4f\n', ...
        roi, N, N_eff, mean(diff_valid), CI95(1), CI95(2), p_two);

    if make_boot_fig
        figure('Name', ['Bootstrap diff: ' roi], 'Position', [120, 120, 700, 450]);
        histogram(diff_valid, 50);
        xline(0, 'k-', 'LineWidth', 2);
        title(sprintf('%s: diff = slope(con) - slope(uncon)', roi));
        xlabel('Bootstrap diff');
        ylabel('Count');
        box off;
    end
end

% BH-FDR across ROIs
if height(boot_summary) > 0
    pvals = boot_summary.p_two(:);
    m = numel(pvals);

    [p_sorted, idx_sort] = sort(pvals, 'ascend');
    ranks = (1:m)';
    q_sorted = p_sorted .* (m ./ ranks);
    q_sorted = flipud(cummin(flipud(q_sorted)));
    q_sorted(q_sorted > 1) = 1;

    qvals = nan(m,1);
    qvals(idx_sort) = q_sorted;

    boot_summary.q_BH    = qvals;
    boot_summary.sig_q05 = double(qvals < 0.05);

    fprintf('\n================ BH-FDR summary across ROIs ================\n');
    disp(boot_summary);
else
    fprintf('\n[Bootstrap] boot_summary is empty: no ROI produced valid results.\n');
end

% Save outputs
% out_boot_mat = fullfile(slope_out_dir, 'bootstrap_slopeDiff_voxelPaired.mat');
% out_boot_csv = fullfile(slope_out_dir, 'bootstrap_slopeDiff_voxelPaired_summary.csv');
% 
% save(out_boot_mat, ...
%     'boot_summary', 'boot_B', 'min_bin_vox', 'rng_seed', ...
%     'boot_slope_uncon_all', 'boot_slope_con_all', 'boot_diff_all', ...
%     'roi_list', 'x_positions_log');
% 
% writetable(boot_summary, out_boot_csv);
% 
% fprintf('\n[Saved] Bootstrap results:\n  MAT: %s\n  CSV: %s\n', out_boot_mat, out_boot_csv);

%% robustness check -- Table S5
% Goal: compare diff_mean and CI95 across different bootstrap sizes and seeds
% This block depends on:
% merged_pref / merged_ve / roi_list / x_positions_log / slope_out_dir

% Sensitivity settings
sens_B_list    = [2000, 5000, 10000];
sens_seed_list = [20260123, 20260124, 20260125];
min_bin_vox    = 2;
make_sens_fig  = false;

% Results table: one row = ROI x B x seed
sens_table = table('Size',[0 9], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'ROI','B','seed','N_vox','N_eff','diff_mean','CI95_low','CI95_high','p_two'});

fprintf('\n================ Sensitivity check starts (diff_mean & CI95) ================\n');
fprintf('B list: %s\n', mat2str(sens_B_list));
fprintf('Seed list: %s\n', mat2str(sens_seed_list));

for r = 1:numel(roi_list)
    roi = roi_list{r};

    if ~isfield(merged_pref, roi) || isempty(merged_pref.(roi).control)
        fprintf('[Sensitivity] ROI=%s has no data. Skipping.\n', roi);
        continue;
    end

    control_data = merged_pref.(roi).control;
    con_data     = merged_pref.(roi).con;
    uncon_data   = merged_pref.(roi).uncon;
    ve_data      = merged_ve.(roi).control;

    valid_idx = (ve_data > 0.3) & ...
                (control_data >= 1.05) & (control_data <= 7.05) & ...
                (con_data     >= 1.05) & (con_data     <= 7.05) & ...
                (uncon_data   >= 1.05) & (uncon_data   <= 7.05);

    control_data = control_data(valid_idx);
    con_data     = con_data(valid_idx);
    uncon_data   = uncon_data(valid_idx);

    N = numel(control_data);
    if N < 10
        fprintf('[Sensitivity] ROI=%s has too few valid voxels (N=%d). Skipping.\n', roi, N);
        continue;
    end

    % Binning according to control PN
    bin_indices = zeros(size(control_data));
    bin_indices(control_data >= 1.05  & control_data <= 1.5) = 1;
    bin_indices(control_data > 1.5 & control_data <= 2.5)   = 2;
    bin_indices(control_data > 2.5 & control_data <= 3.5)   = 3;
    bin_indices(control_data > 3.5 & control_data <= 4.5)   = 4;
    bin_indices(control_data > 4.5 & control_data <= 5.5)   = 5;
    bin_indices(control_data > 5.5 & control_data <= 6.5)   = 6;
    bin_indices(control_data > 6.5 & control_data <= 7.5)   = 7;

    % Precompute log differences
    diff_con   = log(con_data)   - log(control_data);
    diff_uncon = log(uncon_data) - log(control_data);

    for bb = 1:numel(sens_B_list)
        B = sens_B_list(bb);

        for ss = 1:numel(sens_seed_list)
            seed = sens_seed_list(ss);
            rng(seed);

            boot_slope_con   = nan(B,1);
            boot_slope_uncon = nan(B,1);

            for b = 1:B
                boot_idx = randi(N, [N,1]);

                con_bs   = diff_con(boot_idx);
                uncon_bs = diff_uncon(boot_idx);
                bin_bs   = bin_indices(boot_idx);

                con_m   = nan(1,7);
                uncon_m = nan(1,7);

                for bin = 1:7
                    idx = (bin_bs == bin);
                    if sum(idx) < min_bin_vox
                        continue;
                    end
                    con_m(bin)   = mean(con_bs(idx));
                    uncon_m(bin) = mean(uncon_bs(idx));
                end

                vcon   = ~isnan(con_m);
                vuncon = ~isnan(uncon_m);

                if sum(vcon) > 1
                    p = polyfit(x_positions_log(vcon), con_m(vcon), 1);
                    boot_slope_con(b) = p(1);
                end
                if sum(vuncon) > 1
                    p = polyfit(x_positions_log(vuncon), uncon_m(vuncon), 1);
                    boot_slope_uncon(b) = p(1);
                end
            end

            boot_diff = boot_slope_con - boot_slope_uncon;
            valid_iter = ~isnan(boot_diff);
            N_eff = sum(valid_iter);

            if N_eff == 0
                diff_mean = NaN; CI95 = [NaN NaN]; p_two = NaN;
            else
                d = boot_diff(valid_iter);
                diff_mean = mean(d);
                CI95 = prctile(d, [2.5 97.5]);
                p_two = 2 * min(mean(d >= 0), mean(d <= 0));
            end

            sens_table = [sens_table; {string(roi), B, seed, N, N_eff, diff_mean, CI95(1), CI95(2), p_two}]; %#ok<AGROW>

            fprintf('[Sensitivity] ROI=%s | B=%d | seed=%d | N_eff=%d | diff_mean=%.4f | CI95=[%.4f, %.4f]\n', ...
                roi, B, seed, N_eff, diff_mean, CI95(1), CI95(2));
        end
    end
end


roi_names = unique(sens_table.ROI);

stab_table = table('Size',[0 8], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'ROI','n_runs','diff_mean_min','diff_mean_max','diff_mean_range','CI_low_range','CI_high_range','CI_width_range'});

for i = 1:numel(roi_names)
    roi = roi_names(i);
    idx = sens_table.ROI == roi & ~isnan(sens_table.diff_mean);

    if sum(idx) == 0
        continue;
    end

    diff_mean_min = min(sens_table.diff_mean(idx));
    diff_mean_max = max(sens_table.diff_mean(idx));
    diff_mean_rng = diff_mean_max - diff_mean_min;

    CI_low_rng  = max(sens_table.CI95_low(idx))  - min(sens_table.CI95_low(idx));
    CI_high_rng = max(sens_table.CI95_high(idx)) - min(sens_table.CI95_high(idx));

    CI_width = sens_table.CI95_high(idx) - sens_table.CI95_low(idx);
    CI_width_rng = max(CI_width) - min(CI_width);

    stab_table = [stab_table; {roi, sum(idx), diff_mean_min, diff_mean_max, diff_mean_rng, CI_low_rng, CI_high_rng, CI_width_rng}]; %#ok<AGROW>
end

fprintf('\n================ Sensitivity check: per-ROI stability summary (smaller range = more stable) ================\n');
disp(stab_table);

% Optional threshold-based stability classification
diff_thr = 0.02;
ci_thr   = 0.03;

stab_flag = (stab_table.diff_mean_range < diff_thr) & ...
            (stab_table.CI_low_range < ci_thr) & ...
            (stab_table.CI_high_range < ci_thr);

stab_table.stable_by_threshold = double(stab_flag);

fprintf('\n[Decision threshold] diff_mean_range < %.3f and CI endpoint range < %.3f -> stable = 1\n', diff_thr, ci_thr);
disp(stab_table);

% Optional plotting for sensitivity runs
if make_sens_fig
    for i = 1:numel(roi_names)
        roi = roi_names(i);
        idx = sens_table.ROI == roi & ~isnan(sens_table.diff_mean);
        T = sens_table(idx,:);

        if isempty(T); continue; end

        x = 1:height(T);

        figure('Name', ['Sensitivity: ' char(roi)], 'Position', [150, 150, 900, 420]);
        hold on;
        errorbar(x, T.diff_mean, T.diff_mean - T.CI95_low, T.CI95_high - T.diff_mean, 'o', 'LineWidth', 1.2);
        yline(0, 'k-', 'LineWidth', 1.5);
        xlabel('Run index (different B × seed combinations)');
        ylabel('diff_mean and 95% CI');
        title(sprintf('%s: diff_mean and CI95 across different B/seed settings', roi));
        box off;
        hold off;

        fprintf('\n[Plot order] ROI=%s\n', roi);
        disp(T(:,{'B','seed','diff_mean','CI95_low','CI95_high'}));
    end
end

% Save sensitivity results
out_sens_mat = fullfile(slope_out_dir, 'sensitivity_check_diffMean_CI95.mat');
out_sens_csv = fullfile(slope_out_dir, 'sensitivity_check_diffMean_CI95_runs.csv');
out_stab_csv = fullfile(slope_out_dir, 'sensitivity_check_diffMean_CI95_stability.csv');

save(out_sens_mat, 'sens_table', 'stab_table', 'sens_B_list', 'sens_seed_list', 'min_bin_vox');
writetable(sens_table, out_sens_csv);
writetable(stab_table, out_stab_csv);

fprintf('\n[Saved] Sensitivity-check results:\n  MAT: %s\n  CSV (runs): %s\n  CSV (stability): %s\n', ...
    out_sens_mat, out_sens_csv, out_stab_csv);

%% Fig 3 (Panel A+B)
% Goal: based on the log-scale adaptation-effect plot, overlay the bootstrap
% mean fitted line and the 95% CI ribbon.
%
% Requirements in workspace:
% - results.(roi).con_diff_log_mean / sem / uncon_diff_log_mean / sem
% - results.(roi).bootstrap.slope_con / slope_uncon
% - x_positions_log, x_limits, x_labels
% - roi_list
% - preferably intercept_log_6x2 and roi_order
%
% Also generates Fig 3b (Panel 2): point estimate + CI for diff,
% grouped into posterior vs others, with y=0 and BH-FDR stars

% Unified ROI order
roi_order_plot = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};
group_others   = {'NF','NPC1','NPC2','NPC3'};
group_post     = {'NPO','NTO'};

% Colors and plotting parameters
uncon_color = [0.07, 0.62, 1.00];
con_color   = [1.00, 0.00, 0.00];

ribbon_alpha = 0.20;
fit_lw = 2.5;
capw = 0.08;
x_grid = linspace(min(x_positions_log), max(x_positions_log), 200);

% Recompute x_limits / x_labels if missing
if ~exist('x_limits','var') || numel(x_limits)~=2
    x_min = min(x_positions_log); x_max = max(x_positions_log);
    x_range = x_max - x_min; x_margin = 0.05*x_range;
    x_limits = [x_min-x_margin, x_max+x_margin];
end
if ~exist('x_labels','var')
    adapter_numerosity = 40;
    pn_centers = 1:7;
    difference_centers = adapter_numerosity - pn_centers;
    x_labels = arrayfun(@(d) sprintf('%d', d), difference_centers, 'UniformOutput', false);
end

%% Part A: Figure 5A — one bootstrap-fit figure per ROI
for ii = 1:numel(roi_order_plot)
    roi = roi_order_plot{ii};

    if ~isfield(results, roi) || ~isfield(results.(roi),'bootstrap')
        fprintf('[Fig3a] ROI=%s is missing results.(roi).bootstrap. Skipping.\n', roi);
        continue;
    end
    if ~isfield(results.(roi),'con_diff_log_mean') || ~isfield(results.(roi),'uncon_diff_log_mean')
        fprintf('[Fig3a] ROI=%s is missing binned means (con/uncon). Skipping.\n', roi);
        continue;
    end

    % Original binned mean + SEM
    con_m   = results.(roi).con_diff_log_mean;
    uncon_m = results.(roi).uncon_diff_log_mean;
    con_sem   = results.(roi).con_diff_log_sem;
    uncon_sem = results.(roi).uncon_diff_log_sem;

    valid_con_bins   = ~isnan(con_m);
    valid_uncon_bins = ~isnan(uncon_m);

    % Bootstrap slopes
    s_con   = results.(roi).bootstrap.slope_con(:);
    s_uncon = results.(roi).bootstrap.slope_uncon(:);
    valid_it = ~isnan(s_con) & ~isnan(s_uncon);
    s_con   = s_con(valid_it);
    s_uncon = s_uncon(valid_it);

    if numel(s_con) < 50
        fprintf('[Fig3a] ROI=%s has too few valid paired bootstrap iterations (N=%d), but plotting will proceed.\n', roi, numel(s_con));
    end
    if isempty(s_con)
        fprintf('[Fig3a] ROI=%s has no valid bootstrap slopes. Skipping.\n', roi);
        continue;
    end

    % Intercepts: prefer original fit intercepts; otherwise use 0
    b_con = 0; b_uncon = 0;
    if exist('intercept_log_6x2','var') && exist('roi_order','var') && iscell(roi_order)
        ridx = find(strcmp(roi_order, roi), 1);
        if ~isempty(ridx) && size(intercept_log_6x2,2) >= 2
            b_uncon = intercept_log_6x2(ridx,1);
            b_con   = intercept_log_6x2(ridx,2);
        end
    else
        fprintf('[Fig3a] ROI=%s: intercept_log_6x2/roi_order not detected, intercepts are fixed at 0 (CI reflects slope uncertainty only).\n', roi);
    end

    % Generate bootstrap prediction distributions
    Y_con   = bsxfun(@plus, b_con,   bsxfun(@times, s_con,   x_grid));
    Y_uncon = bsxfun(@plus, b_uncon, bsxfun(@times, s_uncon, x_grid));

    y_con_mean   = mean(Y_con,1);
    y_uncon_mean = mean(Y_uncon,1);

    y_con_lo   = prctile(Y_con,   2.5, 1);
    y_con_hi   = prctile(Y_con,  97.5, 1);
    y_uncon_lo = prctile(Y_uncon, 2.5, 1);
    y_uncon_hi = prctile(Y_uncon,97.5, 1);

    % Plot
    figure('Name', ['Fig3a bootstrap fit ROI: ' roi], 'Position', [120, 120, 860, 620]);
    hold on;

    if any(valid_con_bins)
        errorbar(x_positions_log(valid_con_bins), con_m(valid_con_bins), con_sem(valid_con_bins), ...
            'o', 'MarkerSize', 6, 'MarkerFaceColor', con_color, 'MarkerEdgeColor', con_color, ...
            'Color', con_color, 'LineWidth', 1.5, 'CapSize', 8, 'HandleVisibility','off');
    end
    if any(valid_uncon_bins)
        errorbar(x_positions_log(valid_uncon_bins), uncon_m(valid_uncon_bins), uncon_sem(valid_uncon_bins), ...
            'o', 'MarkerSize', 6, 'MarkerFaceColor', uncon_color, 'MarkerEdgeColor', uncon_color, ...
            'Color', uncon_color, 'LineWidth', 1.5, 'CapSize', 8, 'HandleVisibility','off');
    end

    % Bootstrap CI ribbons
    Xp = [x_grid, fliplr(x_grid)];
    Yp = [y_uncon_lo, fliplr(y_uncon_hi)];
    h = patch(Xp, Yp, uncon_color, 'EdgeColor','none', 'HandleVisibility','off');
    try
        set(h,'FaceAlpha',ribbon_alpha);
    catch
        set(h,'FaceColor',0.75*uncon_color + 0.25*[1 1 1]);
    end

    Yp = [y_con_lo, fliplr(y_con_hi)];
    h = patch(Xp, Yp, con_color, 'EdgeColor','none', 'HandleVisibility','off');
    try
        set(h,'FaceAlpha',ribbon_alpha);
    catch
        set(h,'FaceColor',0.75*con_color + 0.25*[1 1 1]);
    end

    % Mean fitted lines
    plot(x_grid, y_uncon_mean, '-', 'Color', uncon_color, 'LineWidth', fit_lw, 'HandleVisibility','off');
    plot(x_grid, y_con_mean,   '-', 'Color', con_color,   'LineWidth', fit_lw, 'HandleVisibility','off');

    % y = 0 reference
    line(x_limits, [0 0], 'Color','k', 'LineWidth', 2, 'HandleVisibility','off');

    ax = gca;
    ax.Box       = 'off';
    ax.TickDir   = 'in';
    ax.LineWidth = 3.0;
    ax.FontSize  = 16;

    xlim(x_limits);
    xticks(x_positions_log);
    xticklabels(x_labels);

    ylim([-1.2, 1.2]);
    yticks(-1.2:0.4:1.2);

    title(['Fig.5A  ' roi '  bootstrap fit (mean \pm 95% CI ribbon)']);
    ylabel('log(PN_{adapt}) - log(PN_{ctrl})');

    h_leg_uncon = plot(nan, nan, 'o', 'MarkerSize', 12, ...
        'MarkerFaceColor', uncon_color, 'MarkerEdgeColor', uncon_color, ...
        'LineStyle','none', 'DisplayName','Unconnected (binned mean \pm SEM)');
    h_leg_con = plot(nan, nan, 'o', 'MarkerSize', 12, ...
        'MarkerFaceColor', con_color, 'MarkerEdgeColor', con_color, ...
        'LineStyle','none', 'DisplayName','Connected (binned mean \pm SEM)');

    legend([h_leg_uncon, h_leg_con], 'Location','best');

    hold off;
end

%% Figure 5B
% Colored light violin per ROI + black mean dot + black 95% CI line
% + fixed ROI order: NF, NPC1-3, NPO, NTO
% + significance stars at fixed y = 0.8
% + narrower violin width
% + fixed y-axis range and ticks: [-0.6, 0.8] with step 0.2

% Sanity checks
if ~exist('boot_diff_all','var') || isempty(boot_diff_all)
    error('boot_diff_all not found or empty. Run the bootstrap section first.');
end
if ~exist('roi_list','var') || isempty(roi_list)
    error('roi_list not found.');
end

% Desired ROI order
desired_order = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};

% Map ROI name to index in boot_diff_all
roi_list_cell = roi_list(:);
idx_map = nan(numel(desired_order),1);
for i = 1:numel(desired_order)
    j = find(strcmp(roi_list_cell, desired_order{i}), 1);
    if ~isempty(j), idx_map(i) = j; end
end

keep = ~isnan(idx_map);
rois_plot = desired_order(keep);
idx_map   = idx_map(keep);

nROI = numel(rois_plot);
if nROI == 0
    error('None of the desired ROIs were found in roi_list: %s', strjoin(desired_order, ', '));
end

% Figure parameters
violin_width    = 0.28;
n_ks_points     = 250;
alpha_violin    = 0.22;
ci_line_width   = 2.6;
mean_marker_sz  = 9;

star_y          = 0.8;
star_fontsize   = 16;

% Fixed y-axis settings
y_lim   = [-0.6, 0.8];
y_ticks = -0.6:0.2:0.8;

% ROI colors
roi_colors = lines(nROI);

% Precompute stats
d_all   = cell(nROI,1);
mu_all  = nan(nROI,1);
CI_all  = nan(nROI,2);
pv_all  = nan(nROI,1);
starstr = repmat({''}, nROI, 1);

for i = 1:nROI
    d = boot_diff_all{idx_map(i)};
    d = d(~isnan(d));
    d_all{i} = d;

    mu_all(i)   = mean(d);
    CI_all(i,:) = prctile(d, [2.5 97.5]);

    p_or_q = NaN;

    if exist('boot_summary','var') && ~isempty(boot_summary) && istable(boot_summary)
        roi_col = boot_summary.ROI;
        if iscell(roi_col)
            match_row = find(strcmp(roi_col, rois_plot{i}), 1);
        else
            match_row = find(roi_col == string(rois_plot{i}), 1);
        end

        if ~isempty(match_row)
            if ismember('q_BH', boot_summary.Properties.VariableNames) && ~isnan(boot_summary.q_BH(match_row))
                p_or_q = boot_summary.q_BH(match_row);
            elseif ismember('p_two', boot_summary.Properties.VariableNames) && ~isnan(boot_summary.p_two(match_row))
                p_or_q = boot_summary.p_two(match_row);
            end
        end
    end

    if isnan(p_or_q)
        p_or_q = 2 * min(mean(d >= 0), mean(d <= 0));
    end

    pv_all(i) = p_or_q;

    if p_or_q < 0.001
        starstr{i} = '***';
    elseif p_or_q < 0.01
        starstr{i} = '**';
    elseif p_or_q < 0.05
        starstr{i} = '*';
    else
        starstr{i} = '';
    end
end

% Plot
figure('Name','Colored violin + mean + 95% CI (diff) | NF NPC1-3 NPO NTO', ...
       'Position',[120 120 1100 460]);
hold on;

yline(0, 'k-', 'LineWidth', 2, 'HandleVisibility','off');

for i = 1:nROI
    d = d_all{i};
    if isempty(d), continue; end

    x0 = i;
    face_col = roi_colors(i,:);

    % Violin density
    use_ks = true;
    try
        [f, yi] = ksdensity(d, 'NumPoints', n_ks_points);
    catch
        use_ks = false;
    end

    if use_ks
        f = f ./ max(f);
        f = f * violin_width;

        x_left  = x0 - f;
        x_right = x0 + f;

        patch([x_left fliplr(x_right)], [yi fliplr(yi)], face_col, ...
            'FaceAlpha', alpha_violin, 'EdgeColor', 'none', 'HandleVisibility','off');

        plot(x_left, yi, '-', 'Color', face_col, 'LineWidth', 0.8, 'HandleVisibility','off');
        plot(x_right, yi, '-', 'Color', face_col, 'LineWidth', 0.8, 'HandleVisibility','off');
    else
        % Fallback if ksdensity is unavailable
        nbins = 70;
        [counts, edges] = histcounts(d, nbins, 'Normalization','pdf');
        centers = edges(1:end-1) + diff(edges)/2;
        counts_sm = smoothdata(counts, 'gaussian', 9);
        counts_sm = counts_sm ./ max(counts_sm);
        counts_sm = counts_sm * violin_width;

        x_left  = x0 - counts_sm;
        x_right = x0 + counts_sm;

        patch([x_left fliplr(x_right)], [centers fliplr(centers)], face_col, ...
            'FaceAlpha', alpha_violin, 'EdgeColor', 'none', 'HandleVisibility','off');

        plot(x_left, centers, '-', 'Color', face_col, 'LineWidth', 0.8, 'HandleVisibility','off');
        plot(x_right, centers, '-', 'Color', face_col, 'LineWidth', 0.8, 'HandleVisibility','off');
    end

    % 95% bootstrap CI
    CI = CI_all(i,:);
    plot([x0 x0], CI, 'k-', 'LineWidth', ci_line_width, 'HandleVisibility','off');

    % Mean marker
    mu = mu_all(i);
    plot(x0, mu, 'o', 'MarkerSize', mean_marker_sz, ...
        'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'HandleVisibility','off');

    % Significance stars
    if ~isempty(starstr{i})
        text(x0, star_y, starstr{i}, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontSize', star_fontsize, 'FontWeight','bold', 'Color','k', ...
            'HandleVisibility','off');
    end
end

% Cosmetics
ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'in';
ax.LineWidth = 3.0;
ax.FontSize  = 16;

xlim([0.5 nROI+0.5]);
xticks(1:nROI);
xticklabels(rois_plot);
xtickangle(30);

ylabel('Bootstrap diff = slope_{con} - slope_{uncon}');
title('ROI-wise condition difference (colored bootstrap violin, mean, 95% CI)');
grid off;

ylim(y_lim);
yticks(y_ticks);

% Minimal legend
h_leg_violin = patch(nan, nan, [0 0 0], 'FaceAlpha', alpha_violin, 'EdgeColor','none', ...
    'DisplayName','Bootstrap distribution (violin)');
h_leg_mean   = plot(nan, nan, 'ko', 'MarkerFaceColor','k', 'MarkerSize', mean_marker_sz, ...
    'DisplayName','Mean');
h_leg_ci     = plot(nan, nan, 'k-', 'LineWidth', ci_line_width, 'DisplayName','95% bootstrap CI');
legend([h_leg_violin, h_leg_mean, h_leg_ci], 'Location','best');

hold off;

%% Figure 6B
% ROI-level three-condition comparison (linear slopes: 40 unconnected + 40 connected + 20 unconnected from Tsouli et al., 2021)
fprintf('\n===== ROI-level (NPC1-3 & NF) three-condition comparison (linear slopes) =====\n');

% Required variables check
if ~exist('slope_lin_6x2','var')
    error('slope_lin_6x2 not found. Please make sure the linear fitting section has been run.');
end

% High numerosity CSV directory
if ~exist('csv_dir','var')
    csv_dir = 'D:\lrs\all-in-one analysis\high adaptation linear';
end

roi_target_3c = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};

% Fixed ROI color mapping
roi_colors = lines(numel(roi_target_3c));
roi_color_map = struct();
for i = 1:numel(roi_target_3c)
    roi_color_map.(roi_target_3c{i}) = roi_colors(i,:);
end

% Read High numerosity slopes (one slope per ROI)
slope_high_map = struct();
for i = 1:numel(roi_target_3c)
    slope_high_map.(roi_target_3c{i}) = NaN;
end

for i = 1:numel(roi_target_3c)
    roiName = roi_target_3c{i};
    csv_file = fullfile(csv_dir, [roiName '.csv']);
    if ~isfile(csv_file)
        warning('Cannot find %s (no High numerosity slope for this ROI).', csv_file);
        continue;
    end

    try
        Tcsv = readtable(csv_file);
        Mcsv = table2array(Tcsv);
    catch
        Mcsv = readmatrix(csv_file);
    end

    if size(Mcsv,2) < 2
        warning('%s has fewer than 2 columns. Skipping.', csv_file);
        continue;
    end

    Mcsv = Mcsv(:,1:2);
    Mcsv = Mcsv(~any(isnan(Mcsv),2),:);
    if size(Mcsv,1) < 7
        warning('%s has fewer than 7 valid rows. Skipping.', csv_file);
        continue;
    end

    yH = Mcsv(1:7,2);
    xH = (1:7)';
    pH = polyfit(xH, yH, 1);
    slope_high_map.(roiName) = pH(1);
end

% Extract three-condition slopes
[tf, idx] = ismember(roi_target_3c, roi_order);
if any(~tf)
    warning('Some entries in roi_target_3c are not in roi_order and will be ignored: %s', strjoin(roi_target_3c(~tf), ', '));
end
idx_sel = idx(tf);
roi_sel = roi_order(idx_sel);

xU_all = slope_lin_6x2(idx_sel, 1);
xC_all = slope_lin_6x2(idx_sel, 2);

xH_all = nan(numel(idx_sel),1);
for i = 1:numel(roi_sel)
    rn = roi_sel{i};
    xH_all(i) = slope_high_map.(rn);
end

mask3 = isfinite(xU_all) & isfinite(xC_all) & isfinite(xH_all);
roi_used_3c = roi_sel(mask3);

xU = xU_all(mask3);
xC = xC_all(mask3);
xH = xH_all(mask3);
N  = numel(xU);

fprintf('Target ROIs: %s\n', strjoin(roi_target_3c, ', '));
fprintf('Number of ROIs with valid slopes in all three conditions: N=%d\n', N);
if N > 0
    fprintf('ROIs used: %s\n', strjoin(roi_used_3c, ', '));
    fprintf('Unconnected slopes: %s\n', mat2str(xU(:)', 5));
    fprintf('Connected slopes: %s\n', mat2str(xC(:)', 5));
    fprintf('High-numerosity slopes: %s\n', mat2str(xH(:)', 5));
else
    error('No ROI has valid slopes in all three conditions. Cannot generate the three-condition paired plot.');
end

% Plot only: no statistical tests
if N >= 1
    cond_names = {'Unconnected','Connected','High numerosity'};
    
    % Increased spacing between bars
    xpos = [1, 2.3, 3.6];
    
    mU = mean(xU); semU = std(xU)/sqrt(N);
    mC = mean(xC); semC = std(xC)/sqrt(N);
    mH = mean(xH); semH = std(xH)/sqrt(N);

    figure('Name','ROI-level slopes: three-condition paired plot (linear slopes)', ...
           'Position',[120,120,850,540]);
    hold on;

    barW = 0.5;
    
    b1 = bar(xpos(1), mU, barW);
    b1.FaceColor = uncon_color;
    b1.EdgeColor = uncon_color;
    
    b2 = bar(xpos(2), mC, barW);
    b2.FaceColor = con_color;
    b2.EdgeColor = con_color;
    
    b3 = bar(xpos(3), mH, barW);
    b3.FaceColor = gray;
    b3.EdgeColor = gray;

    errorbar(xpos(1), mU, semU, 'k', 'LineStyle','none','LineWidth',1.2,'CapSize',10);
    errorbar(xpos(2), mC, semC, 'k', 'LineStyle','none','LineWidth',1.2,'CapSize',10);
    errorbar(xpos(3), mH, semH, 'k', 'LineStyle','none','LineWidth',1.2,'CapSize',10);

    % ROI-specific colors
    for i = 1:N
        rn = roi_used_3c{i};
        if isfield(roi_color_map, rn)
            cROI = roi_color_map.(rn);
        else
            cROI = [0 0 0];
        end

        plot(xpos, [xU(i) xC(i) xH(i)], '-', 'Color',[0.80 0.80 0.80], 'LineWidth',1.2);

        scatter(xpos(1), xU(i), 90, 'o', 'MarkerFaceColor',cROI, 'MarkerEdgeColor','k', 'LineWidth',0.6);
        scatter(xpos(2), xC(i), 90, 'o', 'MarkerFaceColor',cROI, 'MarkerEdgeColor','k', 'LineWidth',0.6);
        scatter(xpos(3), xH(i), 90, 'o', 'MarkerFaceColor',cROI, 'MarkerEdgeColor','k', 'LineWidth',0.6);

        text(xpos(3) + 0.15, xH(i), rn, 'FontSize', 11, 'Interpreter','none');
    end

    yline(0,'k--','LineWidth',1,'HandleVisibility','off');

    set(gca,'XTick',xpos,'XTickLabel',cond_names);
    xlim([xpos(1)-0.7, xpos(3)+0.9]);

    ylabel('Slope (linear-space fit)');
    title(sprintf('ROI-level paired plot, N=%d (no statistical tests)', N));

    ax = gca;
    ax.Box = 'off';
    ax.TickDir = 'in';
    ax.LineWidth = 1.5;
    ax.FontSize = 12;

    % ROI color legend
    hLeg = gobjects(numel(roi_target_3c),1);
    for k = 1:numel(roi_target_3c)
        rn = roi_target_3c{k};
        hLeg(k) = scatter(nan, nan, 90, 'o', 'MarkerFaceColor',roi_color_map.(rn), ...
                          'MarkerEdgeColor','k', 'LineWidth',0.6);
    end
    legend(hLeg, roi_target_3c, 'Location','eastoutside', 'Box','off');

    hold off;
end

% Save three-condition structure
% roi_stats_3cond = struct();
% roi_stats_3cond.roi_target   = roi_target_3c;
% roi_stats_3cond.roi_used     = roi_used_3c;
% roi_stats_3cond.N            = N;
% roi_stats_3cond.slopes_U     = xU;
% roi_stats_3cond.slopes_C     = xC;
% roi_stats_3cond.slopes_H     = xH;
% roi_stats_3cond.note         = 'No statistical tests were performed. ROI colors are fixed per ROI.';

fprintf('ROI-level three-condition plot completed (no statistical tests).\n');

%% Helper functions (end of script)
function q = fdr_bh_simple(p, alpha)
% Input p: column vector; output q: BH-FDR adjusted p values
    if nargin < 2, alpha = 0.05; end %#ok<NASGU>
    p = p(:);
    q = nan(size(p));
    m = sum(isfinite(p));
    if m == 0, return; end

    idxv = find(isfinite(p));
    pv = p(idxv);

    [ps, ord] = sort(pv);
    r = (1:m)';

    adj = ps .* m ./ r;
    adj = flipud(cummin(flipud(adj)));
    adj(adj > 1) = 1;

    qv = nan(size(pv));
    qv(ord) = adj;
    q(idxv) = qv;
end

function stars = p2stars_cn(p)
% Common significance star thresholds for p or q values
    if ~isfinite(p)
        stars = '';
        return;
    end
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = '';
    end
end

function [is_normal, pN, method] = check_normality_simple(x, alpha_norm)
% For small N, normality tests are unstable; this function is only used
% for automatic branching between parametric and nonparametric tests
    x = x(:);
    is_normal = false;
    pN = NaN;
    method = 'none';

    if numel(x) < 3
        method = 'N<3';
        is_normal = false;
        return;
    end

    if exist('lillietest','file') == 2
        try
            [h,p] = lillietest(x, 'Alpha', alpha_norm);
            pN = p; method = 'lillietest';
            is_normal = (h==0);
            return;
        catch
        end
    end

    if exist('jbtest','file') == 2
        try
            [h,p] = jbtest(x, alpha_norm);
            pN = p; method = 'jbtest';
            is_normal = (h==0);
            return;
        catch
        end
    end

    method = 'no_test_available';
    is_normal = false;
end

function v = iqr_simple(x)
% IQR without relying on the Statistics Toolbox: Q3 - Q1
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);
    if n < 2
        v = NaN;
        return;
    end
    x = sort(x);
    q1 = quantile_lin(x, 0.25);
    q3 = quantile_lin(x, 0.75);
    v = q3 - q1;
end

function q = quantile_lin(xs, p)
% xs must be a sorted vector; p in [0,1]
    n = numel(xs);
    if n == 1
        q = xs(1);
        return;
    end
    pos = 1 + (n - 1) * p;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = xs(lo);
    else
        q = xs(lo) + (pos - lo) * (xs(hi) - xs(lo));
    end
end