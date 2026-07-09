%% Figure 4B + TableS1:
%  per-subject slopes + Wilcoxon (signrank) vs 0 + BH-FDR across reported tests
%  + bar plot with FIXED stars at y=0.1
%  + overlay subject dots
%  + ENFORCE PAIRED: if only one condition has slope, drop both
%  + REMOVE OUTLIERS: absolute z-score criterion within each ROI × condition

clear; clc; close all;

%% -------------------- Paths --------------------
% This script assumes the following folder structure:
%
% project_root/
% |-- Figure4B_TableS1_analysis.m
% |-- data/
% |   |-- params/
% |   |   |-- *_merged_data.mat
% |-- results/
%
% If this script is placed in a subfolder such as "analysis", change
% project_root = script_dir;
% to
% project_root = fileparts(script_dir);

script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

project_root = script_dir;

data_dir = fullfile(project_root, 'data', 'params');
out_dir  = fullfile(project_root, 'results');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_mat = fullfile(out_dir, 'slopes_per_subject_merged_Wilcoxon.mat');

% Discover subjects
mat_files = dir(fullfile(data_dir, '*_merged_data.mat'));
if isempty(mat_files)
    error('No *_merged_data.mat files were found in "%s". Please check the path or file naming.', data_dir);
end

[~, ord] = sort({mat_files.name});
mat_files = mat_files(ord);

nSub = numel(mat_files);

subject_id_raw = cell(1, nSub);
subject_list   = cell(1, nSub);

for i = 1:nSub
    [~, name, ~] = fileparts(mat_files(i).name);
    pos = strfind(name, '_merged_data');

    if isempty(pos)
        error('File name "%s" does not match the pattern "*_merged_data.mat".', mat_files(i).name);
    end

    subject_id_raw{i} = name(1:pos-1);
    subject_list{i}   = sprintf('sub%02d', i);
end

fprintf('Found %d subject files in "%s" (anonymous IDs start from sub01):\n', data_dir, nSub);
disp(subject_list);

% Settings
roi_list = {'NF', 'NPC1', 'NPC2', 'NPC3', 'NTO', 'NPO'};
nRoi = numel(roi_list);

adapter_numerosity = 40;
pn_centers = 1:7;
x_bins = log(pn_centers) - log(adapter_numerosity);

VE_TH  = 0.3;
PN_MIN = 1.05;
PN_MAX = 7.05;

MIN_VOX_PER_BIN  = 2;
MIN_BINS_FOR_FIT = 2;

alpha = 0.05;

% Outlier settings
% Absolute z-score criterion:
% within each ROI × condition, remove values satisfying:
% abs((x - mean(x)) / std(x)) > OUTLIER_Z_TH

OUTLIER_Z_TH = 2.5;
MIN_N_FOR_OUTLIER_CHECK = 3;

% Preallocate outputs
slope_con      = nan(nSub, nRoi);
intercept_con  = nan(nSub, nRoi);
r2_con         = nan(nSub, nRoi);

slope_uncon     = nan(nSub, nRoi);
intercept_uncon = nan(nSub, nRoi);
r2_uncon        = nan(nSub, nRoi);

n_vox_used = nan(nSub, nRoi);
for s = 1:nSub
    subj_label = subject_list{s};
    subj_raw   = subject_id_raw{s};

    in_file = fullfile(data_dir, [subj_raw '_merged_data.mat']);

    if ~exist(in_file, 'file')
        warning('Subject file does not exist: %s. Skipping.', in_file);
        continue;
    end

    fprintf('\n=== Processing subject: %s ===\n', subj_label);

    S = load(in_file, 'pref', 've');
    pref = S.pref;
    ve   = S.ve;

    for r = 1:nRoi
        roi = roi_list{r};

        if ~isfield(pref, roi) || ~isfield(ve, roi) || ...
           ~isfield(pref.(roi), 'control') || ~isfield(pref.(roi), 'con') || ...
           ~isfield(pref.(roi), 'uncon')   || ~isfield(ve.(roi), 'control')

            fprintf('  Subject %s is missing ROI %s or related fields. Skipping.\n', subj_label, roi);
            continue;
        end

        control_data = pref.(roi).control(:);
        con_data     = pref.(roi).con(:);
        uncon_data   = pref.(roi).uncon(:);
        ve_data      = ve.(roi).control(:);

        valid_idx = (ve_data > VE_TH) & ...
                    (control_data >= PN_MIN) & (control_data <= PN_MAX) & ...
                    (con_data     >= PN_MIN) & (con_data     <= PN_MAX) & ...
                    (uncon_data   >= PN_MIN) & (uncon_data   <= PN_MAX);

        control_data = control_data(valid_idx);
        con_data     = con_data(valid_idx);
        uncon_data   = uncon_data(valid_idx);

        if numel(control_data) < 5
            fprintf('  ROI %s has too few valid voxels (%d). Skipping.\n', roi, numel(control_data));
            continue;
        end

        fprintf('  ROI %s: number of valid voxels = %d\n', roi, numel(control_data));

        log_control = log(control_data);
        y_con   = log(con_data)   - log_control;
        y_uncon = log(uncon_data) - log_control;

        y_con_bin   = nan(1, 7);
        y_uncon_bin = nan(1, 7);
        n_in_bin    = zeros(1, 7);

        for b = 1:7
            switch b
                case 1
                    idx_bin = (control_data >= 1.05 & control_data <= 1.5);
                case 2
                    idx_bin = (control_data > 1.5 & control_data <= 2.5);
                case 3
                    idx_bin = (control_data > 2.5 & control_data <= 3.5);
                case 4
                    idx_bin = (control_data > 3.5 & control_data <= 4.5);
                case 5
                    idx_bin = (control_data > 4.5 & control_data <= 5.5);
                case 6
                    idx_bin = (control_data > 5.5 & control_data <= 6.5);
                case 7
                    idx_bin = (control_data > 6.5 & control_data <= 7.05);
            end

            n_in_bin(b) = sum(idx_bin);

            if n_in_bin(b) < MIN_VOX_PER_BIN
                continue;
            end

            y_con_bin(b)   = mean(y_con(idx_bin));
            y_uncon_bin(b) = mean(y_uncon(idx_bin));
        end

        n_vox_used(s, r) = sum(n_in_bin);

        % Fit: Connected
        valid_bins_con = (~isnan(y_con_bin)) & (n_in_bin >= MIN_VOX_PER_BIN);

        if sum(valid_bins_con) >= MIN_BINS_FOR_FIT
            x_fit = x_bins(valid_bins_con);
            y_fit = y_con_bin(valid_bins_con);

            p = polyfit(x_fit, y_fit, 1);
            y_hat = polyval(p, x_fit);

            y_mean = mean(y_fit);
            ss_tot = sum((y_fit - y_mean).^2);
            ss_res = sum((y_fit - y_hat).^2);

            r2 = NaN;
            if ss_tot > 0
                r2 = 1 - ss_res / ss_tot;
            end

            slope_con(s, r)     = p(1);
            intercept_con(s, r) = p(2);
            r2_con(s, r)        = r2;
        end

        % Fit: Unconnected
        valid_bins_uncon = (~isnan(y_uncon_bin)) & (n_in_bin >= MIN_VOX_PER_BIN);

        if sum(valid_bins_uncon) >= MIN_BINS_FOR_FIT
            x_fit = x_bins(valid_bins_uncon);
            y_fit = y_uncon_bin(valid_bins_uncon);

            p = polyfit(x_fit, y_fit, 1);
            y_hat = polyval(p, x_fit);

            y_mean = mean(y_fit);
            ss_tot = sum((y_fit - y_mean).^2);
            ss_res = sum((y_fit - y_hat).^2);

            r2 = NaN;
            if ss_tot > 0
                r2 = 1 - ss_res / ss_tot;
            end

            slope_uncon(s, r)     = p(1);
            intercept_uncon(s, r) = p(2);
            r2_uncon(s, r)        = r2;
        end
    end
end

% Save raw slopes
slope_con_raw       = slope_con;
slope_uncon_raw     = slope_uncon;
intercept_con_raw   = intercept_con;
intercept_uncon_raw = intercept_uncon;
r2_con_raw          = r2_con;
r2_uncon_raw        = r2_uncon;

% Remove outliers by absolute z-score threshold
fprintf('\n===== Remove outliers using mean ± %.1f SD within each ROI × condition =====\n', OUTLIER_Z_TH);

outlier_log = struct( ...
    'subject', {}, ...
    'roi', {}, ...
    'condition', {}, ...
    'value', {}, ...
    'mean_before', {}, ...
    'sd_before', {}, ...
    'abs_z', {} );

for r = 1:nRoi
    roi = roi_list{r};

    % Connected condition
    x = slope_con(:, r);
    valid = isfinite(x);
    n_valid = sum(valid);

    if n_valid >= MIN_N_FOR_OUTLIER_CHECK
        mu = mean(x(valid));
        sd = std(x(valid));

        if isfinite(sd) && sd > 0
            z_abs = abs((x - mu) ./ sd);
            outlier_mask_con = valid & (z_abs > OUTLIER_Z_TH);

            if any(outlier_mask_con)
                idx = find(outlier_mask_con);

                for k = 1:numel(idx)
                    sidx = idx(k);

                    fprintf('Remove outlier: subject=%s, ROI=%s, condition=Con, slope=%.4f, mean=%.4f, SD=%.4f, |z|=%.4f\n', ...
                        subject_list{sidx}, roi, slope_con(sidx, r), mu, sd, z_abs(sidx));

                    outlier_log(end+1).subject     = subject_list{sidx};
                    outlier_log(end).roi           = roi;
                    outlier_log(end).condition     = 'Con';
                    outlier_log(end).value         = slope_con(sidx, r);
                    outlier_log(end).mean_before   = mu;
                    outlier_log(end).sd_before     = sd;
                    outlier_log(end).abs_z         = z_abs(sidx);
                end

                slope_con(outlier_mask_con, r)     = NaN;
                intercept_con(outlier_mask_con, r) = NaN;
                r2_con(outlier_mask_con, r)        = NaN;
            else
                fprintf('ROI %s, Con: no outlier detected. N=%d, mean=%.4f, SD=%.4f\n', ...
                    roi, n_valid, mu, sd);
            end
        else
            fprintf('ROI %s, Con: SD is zero or invalid. Skip outlier detection.\n', roi);
        end
    else
        fprintf('ROI %s, Con: too few valid subjects for outlier detection. N=%d\n', roi, n_valid);
    end

    % Unconnected condition
    x = slope_uncon(:, r);
    valid = isfinite(x);
    n_valid = sum(valid);

    if n_valid >= MIN_N_FOR_OUTLIER_CHECK
        mu = mean(x(valid));
        sd = std(x(valid));

        if isfinite(sd) && sd > 0
            z_abs = abs((x - mu) ./ sd);
            outlier_mask_uncon = valid & (z_abs > OUTLIER_Z_TH);

            if any(outlier_mask_uncon)
                idx = find(outlier_mask_uncon);

                for k = 1:numel(idx)
                    sidx = idx(k);

                    fprintf('Remove outlier: subject=%s, ROI=%s, condition=Uncon, slope=%.4f, mean=%.4f, SD=%.4f, |z|=%.4f\n', ...
                        subject_list{sidx}, roi, slope_uncon(sidx, r), mu, sd, z_abs(sidx));

                    outlier_log(end+1).subject     = subject_list{sidx};
                    outlier_log(end).roi           = roi;
                    outlier_log(end).condition     = 'Uncon';
                    outlier_log(end).value         = slope_uncon(sidx, r);
                    outlier_log(end).mean_before   = mu;
                    outlier_log(end).sd_before     = sd;
                    outlier_log(end).abs_z         = z_abs(sidx);
                end

                slope_uncon(outlier_mask_uncon, r)     = NaN;
                intercept_uncon(outlier_mask_uncon, r) = NaN;
                r2_uncon(outlier_mask_uncon, r)        = NaN;
            else
                fprintf('ROI %s, Uncon: no outlier detected. N=%d, mean=%.4f, SD=%.4f\n', ...
                    roi, n_valid, mu, sd);
            end
        else
            fprintf('ROI %s, Uncon: SD is zero or invalid. Skip outlier detection.\n', roi);
        end
    else
        fprintf('ROI %s, Uncon: too few valid subjects for outlier detection. N=%d\n', roi, n_valid);
    end
end

fprintf('\nTotal removed outliers before paired enforcement: %d\n', numel(outlier_log));

if ~isempty(outlier_log)
    outlier_table = struct2table(outlier_log);
    disp(outlier_table);
else
    outlier_table = table();
end

% ENFORCE PAIRED
fprintf('\n===== Enforce paired slopes per ROI after outlier removal =====\n');

for r = 1:nRoi
    roi = roi_list{r};

    has_con   = isfinite(slope_con(:, r));
    has_uncon = isfinite(slope_uncon(:, r));

    drop_mask = xor(has_con, has_uncon);

    if any(drop_mask)
        fprintf('ROI %s: %d subjects had only one condition after outlier removal; both-condition data are removed:\n', ...
            roi, sum(drop_mask));

        disp(subject_list(drop_mask));

        slope_con(drop_mask, r)      = NaN;
        intercept_con(drop_mask, r)  = NaN;
        r2_con(drop_mask, r)         = NaN;

        slope_uncon(drop_mask, r)     = NaN;
        intercept_uncon(drop_mask, r) = NaN;
        r2_uncon(drop_mask, r)        = NaN;

        n_vox_used(drop_mask, r)      = NaN;
    else
        fprintf('ROI %s: no paired removal needed.\n', roi);
    end
end

%% Statistics:
%  (A) Con vs 0
%  (B) Uncon vs 0
%  + BH-FDR across reported slope-vs-zero tests only

fprintf('\n===== Wilcoxon signrank after outlier removal =====\n');

p_con   = nan(1, nRoi);
p_uncon = nan(1, nRoi);
z_con   = nan(1, nRoi);
z_uncon = nan(1, nRoi);
n_con   = nan(1, nRoi);
n_uncon = nan(1, nRoi);

median_con    = nan(1, nRoi);
median_uncon  = nan(1, nRoi);
iqr_con       = nan(1, nRoi);
iqr_uncon     = nan(1, nRoi);

for r = 1:nRoi
    roi = roi_list{r};

    % ---------- Con vs 0 ----------
    x = slope_con(:, r);
    x = x(isfinite(x));

    n_con(r) = numel(x);

    if n_con(r) >= 2
        median_con(r) = median(x);
        iqr_con(r)    = iqr_simple(x);

        try
            [p, ~, stats] = signrank(x, 0, 'method', 'exact');
        catch
            [p, ~, stats] = signrank(x, 0);
        end

        p_con(r) = p;

        if isstruct(stats) && isfield(stats, 'zval')
            z_con(r) = stats.zval;
        end

        fprintf('ROI %s, Con vs 0:      N=%d, median=%.4f, IQR=%.4f, p=%.4g\n', ...
            roi, n_con(r), median_con(r), iqr_con(r), p_con(r));
    else
        fprintf('ROI %s, Con vs 0:      N<2, skip. N=%d\n', roi, n_con(r));
    end

    % ---------- Uncon vs 0 ----------
    x = slope_uncon(:, r);
    x = x(isfinite(x));

    n_uncon(r) = numel(x);

    if n_uncon(r) >= 2
        median_uncon(r) = median(x);
        iqr_uncon(r)    = iqr_simple(x);

        try
            [p, ~, stats] = signrank(x, 0, 'method', 'exact');
        catch
            [p, ~, stats] = signrank(x, 0);
        end

        p_uncon(r) = p;

        if isstruct(stats) && isfield(stats, 'zval')
            z_uncon(r) = stats.zval;
        end

        fprintf('ROI %s, Uncon vs 0:    N=%d, median=%.4f, IQR=%.4f, p=%.4g\n', ...
            roi, n_uncon(r), median_uncon(r), iqr_uncon(r), p_uncon(r));
    else
        fprintf('ROI %s, Uncon vs 0:    N<2, skip. N=%d\n', roi, n_uncon(r));
    end
end

% BH-FDR across reported slope-vs-zero tests
all_p = [p_con(:); p_uncon(:)];
valid_p = isfinite(all_p);

[~, crit_p, adj_p_valid] = simple_fdr_bh(all_p(valid_p), alpha);

adj_p_full = nan(size(all_p));
adj_p_full(valid_p) = adj_p_valid;

p_con_fdr   = adj_p_full(1:nRoi);
p_uncon_fdr = adj_p_full(nRoi+1:2*nRoi);

sig_con_fdr   = isfinite(p_con_fdr)   & (p_con_fdr   < alpha);
sig_uncon_fdr = isfinite(p_uncon_fdr) & (p_uncon_fdr < alpha);

fprintf('\n===== BH-FDR across reported slope-vs-zero tests: alpha=%.3f =====\n', alpha);
fprintf('crit_p = %.4g\n', crit_p);

for r = 1:nRoi
    if isfinite(p_con(r))
        fprintf('ROI %s Con vs0:     N=%d, median=%.4f, IQR=%.4f, p=%.4g, q=%.4g, sig=%d\n', ...
            roi_list{r}, n_con(r), median_con(r), iqr_con(r), p_con(r), p_con_fdr(r), sig_con_fdr(r));
    end

    if isfinite(p_uncon(r))
        fprintf('ROI %s Uncon vs0:   N=%d, median=%.4f, IQR=%.4f, p=%.4g, q=%.4g, sig=%d\n', ...
            roi_list{r}, n_uncon(r), median_uncon(r), iqr_uncon(r), p_uncon(r), p_uncon_fdr(r), sig_uncon_fdr(r));
    end
end

%% -------------------- Save --------------------
% save(out_mat, ...
%     'subject_list', 'subject_id_raw', 'roi_list', ...
%     'slope_con', 'slope_uncon', ...
%     'intercept_con', 'intercept_uncon', ...
%     'r2_con', 'r2_uncon', ...
%     'slope_con_raw', 'slope_uncon_raw', ...
%     'intercept_con_raw', 'intercept_uncon_raw', ...
%     'r2_con_raw', 'r2_uncon_raw', ...
%     'n_vox_used', ...
%     'pn_centers', 'x_bins', 'adapter_numerosity', ...
%     'VE_TH', 'PN_MIN', 'PN_MAX', ...
%     'MIN_VOX_PER_BIN', 'MIN_BINS_FOR_FIT', ...
%     'OUTLIER_Z_TH', 'MIN_N_FOR_OUTLIER_CHECK', ...
%     'outlier_log', 'outlier_table', ...
%     'alpha', ...
%     'p_con', 'p_uncon', ...
%     'p_con_fdr', 'p_uncon_fdr', ...
%     'sig_con_fdr', 'sig_uncon_fdr', 'crit_p', ...
%     'z_con', 'z_uncon', ...
%     'n_con', 'n_uncon', ...
%     'median_con', 'median_uncon', ...
%     'iqr_con', 'iqr_uncon');
% 
% fprintf('\nAll processing completed. Results saved to:\n%s\n', out_mat);

% Plot
mean_con   = nan(1, nRoi);
sem_con    = nan(1, nRoi);
mean_uncon = nan(1, nRoi);
sem_uncon  = nan(1, nRoi);

for r = 1:nRoi
    sl_con   = slope_con(:, r);
    sl_uncon = slope_uncon(:, r);

    valid_con   = isfinite(sl_con);
    valid_uncon = isfinite(sl_uncon);

    if any(valid_con)
        mean_con(r) = mean(sl_con(valid_con));
        sem_con(r)  = std(sl_con(valid_con)) / sqrt(sum(valid_con));
    end

    if any(valid_uncon)
        mean_uncon(r) = mean(sl_uncon(valid_uncon));
        sem_uncon(r)  = std(sl_uncon(valid_uncon)) / sqrt(sum(valid_uncon));
    end
end

figure('Name', 'Figure 4B', ...
       'Position', [100, 100, 1100, 520]);

cond_colors = [0.07, 0.62, 1.00; ...
               1.00, 0.00, 0.00];

uncon_color = cond_colors(1, :);
con_color   = cond_colors(2, :);

x0    = 1:nRoi;
delta = 0.18;
barW  = 0.28;

x_con   = x0 - delta;
x_uncon = x0 + delta;

hold on;

bh = gobjects(1, 2);

bh(1) = bar(x_con(:),   mean_con(:),   barW);
bh(2) = bar(x_uncon(:), mean_uncon(:), barW);

bh(1).DisplayName = 'Connected (Con)';
bh(2).DisplayName = 'Unconnected (Uncon)';

bh(1).FaceColor = con_color;
bh(2).FaceColor = uncon_color;

bh(1).EdgeColor = con_color;
bh(2).EdgeColor = uncon_color;

xtips_con   = x_con;
xtips_uncon = x_uncon;

hE1 = errorbar(xtips_con, mean_con(:), sem_con(:), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 8);

hE2 = errorbar(xtips_uncon, mean_uncon(:), sem_uncon(:), ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 8);

set([hE1 hE2], 'HandleVisibility', 'off');

h0 = yline(0, 'k-', 'LineWidth', 2);
set(h0, 'HandleVisibility', 'off');

set(gca, 'XTick', x0, 'XTickLabel', roi_list, 'XTickLabelRotation', 0);

xlabel('Maps');
ylabel('Mean Slope Values');
title('Figure 4B');

ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'in';
ax.LineWidth = 1.5;
ax.FontSize  = 12;
grid off;

xlim([0.5, nRoi + 0.5]);
ylim([-1.2, 0.5]);

% Overlay subject dots
dot_color = [0 0 0];
dot_size = 18;

if nSub == 1
    jitter = 0;
else
    jitter_span = 0.18;
    jitter = linspace(-jitter_span / 2, jitter_span / 2, nSub);
end

for s = 1:nSub
    y = slope_con(s, :);
    m = isfinite(y);

    if any(m)
        scatter(xtips_con(m) + jitter(s), y(m), dot_size, 'o', ...
            'MarkerFaceColor', dot_color, ...
            'MarkerEdgeColor', dot_color, ...
            'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    y = slope_uncon(s, :);
    m = isfinite(y);

    if any(m)
        scatter(xtips_uncon(m) + jitter(s), y(m), dot_size, 'o', ...
            'MarkerFaceColor', dot_color, ...
            'MarkerEdgeColor', dot_color, ...
            'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end
end

legend([bh(1), bh(2)], {'Connected (Con)', 'Unconnected (Uncon)'}, ...
    'Location', 'bestoutside', 'Interpreter', 'none');
fprintf('\n=== FDR-corrected q values and star labels after outlier removal ===\n');

y_star_fixed = 0.1;
yL = ylim;
yR = yL(2) - yL(1);

y_sig_line = y_star_fixed;
y_sig_text = y_sig_line + 0.03 * yR;

if y_sig_text > yL(2) - 0.01 * yR
    y_sig_text = yL(2) - 0.01 * yR;
end

halfLineW = barW / 2;

for r = 1:nRoi
    % Con vs 0
    if isfinite(p_con_fdr(r))
        stars_con = p2stars(p_con_fdr(r));

        fprintf('ROI %s, Con vs0:     q = %.4g, stars = "%s"\n', ...
            roi_list{r}, p_con_fdr(r), stars_con);

        if ~isempty(stars_con)
            hL = line([xtips_con(r) - halfLineW, xtips_con(r) + halfLineW], ...
                      [y_sig_line, y_sig_line], ...
                      'Color', [0 0 0], ...
                      'LineWidth', 1.5, ...
                      'Clipping', 'on');

            set(hL, 'HandleVisibility', 'off');

            text(xtips_con(r), y_sig_text, stars_con, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'Color', [0 0 0], ...
                'FontSize', 18, ...
                'FontWeight', 'bold', ...
                'Clipping', 'on');
        end
    end

    % Uncon vs 0
    if isfinite(p_uncon_fdr(r))
        stars_uncon = p2stars(p_uncon_fdr(r));

        fprintf('ROI %s, Uncon vs0:   q = %.4g, stars = "%s"\n', ...
            roi_list{r}, p_uncon_fdr(r), stars_uncon);

        if ~isempty(stars_uncon)
            hL = line([xtips_uncon(r) - halfLineW, xtips_uncon(r) + halfLineW], ...
                      [y_sig_line, y_sig_line], ...
                      'Color', [0 0 0], ...
                      'LineWidth', 1.5, ...
                      'Clipping', 'on');

            set(hL, 'HandleVisibility', 'off');

            text(xtips_uncon(r), y_sig_text, stars_uncon, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'Color', [0 0 0], ...
                'FontSize', 18, ...
                'FontWeight', 'bold', ...
                'Clipping', 'on');
        end
    end
end

hold off;

%% Helper functions
function stars = p2stars(p)
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

function [h, crit_p, adj_p] = simple_fdr_bh(pvals, alpha)
    if nargin < 2
        alpha = 0.05;
    end

    p = pvals(:);
    m = numel(p);

    if m == 0
        h = [];
        crit_p = NaN;
        adj_p = [];
        return;
    end

    [ps, idx] = sort(p);
    r = (1:m)';

    thresh = r / m * alpha;
    below  = find(ps <= thresh);

    h = false(m, 1);

    if ~isempty(below)
        k = below(end);
        crit_p = ps(k);
        h(1:k) = true;
    else
        crit_p = 0;
    end

    adj_ps = ps .* m ./ r;
    adj_ps = flipud(cummin(flipud(adj_ps)));
    adj_ps(adj_ps > 1) = 1;

    adj_p = nan(m, 1);
    adj_p(idx) = adj_ps;
end

function v = iqr_simple(x)
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