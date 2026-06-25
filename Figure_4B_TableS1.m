%% Figure 4B + TableS1:
%  per-subject slopes + Wilcoxon (signrank) vs 0 + BH-FDR
%  + bar plot with FIXED stars at y=0.1
%  + overlay subject dots (NOW: all black)
%  + ENFORCE PAIRED: if only one condition has slope, drop both (set NaN)
%  + ADD paired test: Con vs Uncon (signrank) + BH-FDR
clear; clc; close all;

%% -------------------- Paths --------------------
data_dir = 'D:\lrs\all-in-one analysis\params';
out_mat  = fullfile(data_dir, 'slopes_per_subject_merged_Wilcoxon.mat');

%% -------------------- Discover subjects --------------------
mat_files = dir(fullfile(data_dir, '*_merged_data.mat'));
if isempty(mat_files)
    error('No *_merged_data.mat files were found in "%s". Please check the path or file naming.', data_dir);
end

% Sort by file name to ensure a consistent subject order across runs
[~, ord] = sort({mat_files.name});
mat_files = mat_files(ord);

nSub = numel(mat_files);

% Real file prefix (from file name), used for loading files
subject_id_raw = cell(1, nSub);

% Anonymous subject IDs for display/saving/statistics: starting from sub01
subject_list = cell(1, nSub);

for i = 1:nSub
    [~, name, ~] = fileparts(mat_files(i).name);
    pos = strfind(name, '_merged_data');
    if isempty(pos)
        error('File name "%s" does not match the pattern "*_merged_data.mat".', mat_files(i).name);
    end
    subject_id_raw{i} = name(1:pos-1);               % used for file loading
    subject_list{i}   = sprintf('sub%02d', i);       % used for display/saving (starting from sub01)
end

fprintf('Found %d subject files in "%s" (anonymous IDs start from sub01):\n', data_dir, nSub);
disp(subject_list);

%% -------------------- Settings --------------------
roi_list = {'NF', 'NPC1', 'NPC2', 'NPC3', 'NTO', 'NPO' };
nRoi = numel(roi_list);

adapter_numerosity = 40;
pn_centers = 1:7;
x_bins = log(pn_centers) - log(adapter_numerosity);  % x = log(PN_control_bin_center) - log(adapter)

% Filter
VE_TH  = 0.3;
PN_MIN = 1.05;
PN_MAX = 7.05;

MIN_VOX_PER_BIN  = 2;
MIN_BINS_FOR_FIT = 2;

alpha = 0.05;

%% -------------------- Preallocate outputs --------------------
slope_con      = nan(nSub, nRoi);
intercept_con  = nan(nSub, nRoi);
r2_con         = nan(nSub, nRoi);

slope_uncon     = nan(nSub, nRoi);
intercept_uncon = nan(nSub, nRoi);
r2_uncon        = nan(nSub, nRoi);

n_vox_used = nan(nSub, nRoi);

%% -------------------- Main loop: subject × ROI (MERGED voxels) --------------------
for s = 1:nSub
    subj_label = subject_list{s};   % sub01, sub02, ...
    subj_raw   = subject_id_raw{s}; % actual file prefix

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

        % ---------- bin by PN_control ----------
        y_con_bin   = nan(1, 7);
        y_uncon_bin = nan(1, 7);
        n_in_bin    = zeros(1, 7);

        for b = 1:7
            switch b
                case 1
                    idx_bin = (control_data >= 1.05   & control_data <= 1.5);
                case 2
                    idx_bin = (control_data >  1.5 & control_data <= 2.5);
                case 3
                    idx_bin = (control_data >  2.5 & control_data <= 3.5);
                case 4
                    idx_bin = (control_data >  3.5 & control_data <= 4.5);
                case 5
                    idx_bin = (control_data >  4.5 & control_data <= 5.5);
                case 6
                    idx_bin = (control_data >  5.5 & control_data <= 6.5);
                case 7
                    idx_bin = (control_data >  6.5 & control_data <= 7.05);
            end

            n_in_bin(b) = sum(idx_bin);
            if n_in_bin(b) < MIN_VOX_PER_BIN
                continue;
            end
            y_con_bin(b)   = mean(y_con(idx_bin));
            y_uncon_bin(b) = mean(y_uncon(idx_bin));
        end

        n_vox_used(s, r) = sum(n_in_bin);

        % ---------- Fit: Con ----------
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
                r2 = 1 - ss_res/ss_tot;
            end

            slope_con(s, r)     = p(1);
            intercept_con(s, r) = p(2);
            r2_con(s, r)        = r2;
        end

        % ---------- Fit: Uncon ----------
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
                r2 = 1 - ss_res/ss_tot;
            end

            slope_uncon(s, r)     = p(1);
            intercept_uncon(s, r) = p(2);
            r2_uncon(s, r)        = r2;
        end
    end
end

%%  ENFORCE PAIRED
fprintf('\n===== Enforce paired slopes per ROI (Con & Uncon must both exist) =====\n');
for r = 1:nRoi
    roi = roi_list{r};

    has_con   = isfinite(slope_con(:, r));
    has_uncon = isfinite(slope_uncon(:, r));

    drop_mask = xor(has_con, has_uncon);

    if any(drop_mask)
        fprintf('ROI %s: %d subjects had only one condition fitted; both-condition data for this ROI have been forcibly removed:\n', roi, sum(drop_mask));
        disp(subject_list(drop_mask));

        slope_con(drop_mask, r)      = NaN;
        intercept_con(drop_mask, r)  = NaN;
        r2_con(drop_mask, r)         = NaN;

        slope_uncon(drop_mask, r)     = NaN;
        intercept_uncon(drop_mask, r) = NaN;
        r2_uncon(drop_mask, r)        = NaN;

        n_vox_used(drop_mask, r)      = NaN;
    else
        fprintf('ROI %s: no removal needed (valid subjects are fully paired for Con/Uncon).\n', roi);
    end
end

%% =====================================================================
%  (A) Single condition vs 0: Wilcoxon signrank
%  (B) Condition difference: paired Wilcoxon signrank (Con vs Uncon)
%  + BH-FDR: corrected across 3*nRoi tests together
%  + Print: median, IQR, N, raw p, q(FDR)

fprintf('\n===== Wilcoxon signrank: (1) Con vs 0, (2) Uncon vs 0, (3) Con vs Uncon =====\n');

% --- (A) vs 0 ---
p_con   = nan(1, nRoi);
p_uncon = nan(1, nRoi);
z_con   = nan(1, nRoi);
z_uncon = nan(1, nRoi);
n_con   = nan(1, nRoi);
n_uncon = nan(1, nRoi);

median_con  = nan(1, nRoi);
median_uncon= nan(1, nRoi);
iqr_con     = nan(1, nRoi);
iqr_uncon   = nan(1, nRoi);

% --- (B) paired difference ---
p_diff  = nan(1, nRoi);   % Con vs Uncon
z_diff  = nan(1, nRoi);
n_diff  = nan(1, nRoi);
median_diff = nan(1, nRoi);
iqr_diff    = nan(1, nRoi);

for r = 1:nRoi
    roi = roi_list{r};

    % Con vs 0
    x = slope_con(:, r);
    x = x(isfinite(x));
    n_con(r) = numel(x);
    if n_con(r) >= 2
        median_con(r) = median(x);
        iqr_con(r)    = iqr_simple(x);

        try
            [p,~,stats] = signrank(x, 0, 'method', 'exact');
        catch
            [p,~,stats] = signrank(x, 0);
        end
        p_con(r) = p;
        if isstruct(stats) && isfield(stats,'zval'), z_con(r) = stats.zval; end

        fprintf('ROI %s, Con vs 0:      N=%d, median=%.4f, IQR=%.4f, p=%.4g\n', ...
            roi, n_con(r), median_con(r), iqr_con(r), p_con(r));
    else
        fprintf('ROI %s, Con vs 0:      N<2 (N=%d), skip\n', roi, n_con(r));
    end

    % Uncon vs 0
    x = slope_uncon(:, r);
    x = x(isfinite(x));
    n_uncon(r) = numel(x);
    if n_uncon(r) >= 2
        median_uncon(r) = median(x);
        iqr_uncon(r)    = iqr_simple(x);

        try
            [p,~,stats] = signrank(x, 0, 'method', 'exact');
        catch
            [p,~,stats] = signrank(x, 0);
        end
        p_uncon(r) = p;
        if isstruct(stats) && isfield(stats,'zval'), z_uncon(r) = stats.zval; end

        fprintf('ROI %s, Uncon vs 0:    N=%d, median=%.4f, IQR=%.4f, p=%.4g\n', ...
            roi, n_uncon(r), median_uncon(r), iqr_uncon(r), p_uncon(r));
    else
        fprintf('ROI %s, Uncon vs 0:    N<2 (N=%d), skip\n', roi, n_uncon(r));
    end

    % Con vs Uncon (paired signrank) —— perform signrank(d,0) on the difference d
    xc = slope_con(:, r);
    xu = slope_uncon(:, r);
    pair_mask = isfinite(xc) & isfinite(xu);   % enforce pairing
    d = xc(pair_mask) - xu(pair_mask);
    n_diff(r) = numel(d);
    if n_diff(r) >= 2
        median_diff(r) = median(d);
        iqr_diff(r)    = iqr_simple(d);

        try
            [p,~,stats] = signrank(d, 0, 'method', 'exact');
        catch
            [p,~,stats] = signrank(d, 0);
        end
        p_diff(r) = p;
        if isstruct(stats) && isfield(stats,'zval'), z_diff(r) = stats.zval; end

        fprintf('ROI %s, Con-Uncon:     N=%d, median(d)=%.4f, IQR(d)=%.4f, p=%.4g\n', ...
            roi, n_diff(r), median_diff(r), iqr_diff(r), p_diff(r));
    else
        fprintf('ROI %s, Con-Uncon:     N<2 (N=%d), skip\n', roi, n_diff(r));
    end
end

%% -------------------- BH-FDR across 3*nRoi tests --------------------
all_p = [p_con(:); p_uncon(:); p_diff(:)];
valid_p = isfinite(all_p);

[~, crit_p, adj_p_valid] = simple_fdr_bh(all_p(valid_p), alpha);

adj_p_full = nan(size(all_p));
adj_p_full(valid_p) = adj_p_valid;

% Here p_*_fdr are actually the FDR-corrected q values (adjusted p values)
p_con_fdr   = adj_p_full(1:nRoi);
p_uncon_fdr = adj_p_full(nRoi+1:2*nRoi);
p_diff_fdr  = adj_p_full(2*nRoi+1:3*nRoi);

% Define significance using q<alpha (to avoid inconsistency between sig and q caused by back-filling h)
sig_con_fdr   = isfinite(p_con_fdr)   & (p_con_fdr   < alpha);
sig_uncon_fdr = isfinite(p_uncon_fdr) & (p_uncon_fdr < alpha);
sig_diff_fdr  = isfinite(p_diff_fdr)  & (p_diff_fdr  < alpha);

fprintf('\n===== BH-FDR across (Con vs0, Uncon vs0, Con-Uncon) : alpha=%.3f =====\n', alpha);
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
    if isfinite(p_diff(r))
        fprintf('ROI %s Con-Uncon:   N=%d, median(d)=%.4f, IQR(d)=%.4f, p=%.4g, q=%.4g, sig=%d\n', ...
            roi_list{r}, n_diff(r), median_diff(r), iqr_diff(r), p_diff(r), p_diff_fdr(r), sig_diff_fdr(r));
    end
end

%% -------------------- Save --------------------
save(out_mat, ...
    'subject_list','roi_list', ...
    'slope_con','slope_uncon', ...
    'intercept_con','intercept_uncon', ...
    'r2_con','r2_uncon', ...
    'n_vox_used', ...
    'pn_centers','x_bins','adapter_numerosity', ...
    'alpha', ...
    'p_con','p_uncon','p_diff', ...
    'p_con_fdr','p_uncon_fdr','p_diff_fdr', ...
    'sig_con_fdr','sig_uncon_fdr','sig_diff_fdr','crit_p', ...
    'z_con','z_uncon','z_diff', ...
    'n_con','n_uncon','n_diff', ...
    'median_con','median_uncon','median_diff', ...
    'iqr_con','iqr_uncon','iqr_diff');

fprintf('\nAll processing completed. Results have been saved to: %s\n', out_mat);

%% Plot (manual x offsets + fixed ylim + sig lines)
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

figure('Name', 'MERGED-voxel slopes with FDR stars (Wilcoxon, fixed y=0.1)', 'Position', [100, 100, 1100, 520]);

% Colors: row 1 = Unconnected (blue), row 2 = Connected (red)
cond_colors = [0.07, 0.62, 1.00;   % Unconnected -> blue
               1.00, 0.00, 0.00];  % Connected   -> red
uncon_color = cond_colors(1,:);    % Unconnected
con_color   = cond_colors(2,:);    % Connected

% Manually assign x offsets for the two conditions within each ROI so they are slightly separated
x0    = 1:nRoi;     % center position of each ROI (spacing = 1)
delta = 0.18;       % left-right offset for the two conditions (2*delta < 1)
barW  = 0.28;       % bar width (must be < 2*delta to avoid overlap)

x_con   = x0 - delta;   % x positions for Connected
x_uncon = x0 + delta;   % x positions for Unconnected

hold on;

% Draw the two bar groups separately (to better control within-ROI spacing)
bh = gobjects(1,2);
bh(1) = bar(x_con(:),   mean_con(:),   barW);
bh(2) = bar(x_uncon(:), mean_uncon(:), barW);

bh(1).DisplayName = 'Connected (Con)';
bh(2).DisplayName = 'Unconnected (Uncon)';

% Bar colors
bh(1).FaceColor = con_color;
bh(2).FaceColor = uncon_color;
bh(1).EdgeColor = con_color;
bh(2).EdgeColor = uncon_color;

% Error bar positions (aligned with bar centers)
xtips_con   = x_con;
xtips_uncon = x_uncon;

hE1 = errorbar(xtips_con,   mean_con(:),   sem_con(:),   'k', 'LineStyle','none','LineWidth',1.0,'CapSize',8);
hE2 = errorbar(xtips_uncon, mean_uncon(:), sem_uncon(:), 'k', 'LineStyle','none','LineWidth',1.0,'CapSize',8);
set([hE1 hE2], 'HandleVisibility','off');  % prevent data1/data2 from appearing in the legend

% y=0 reference line: solid black and thicker (LineWidth=2), not shown in legend
h0 = yline(0, 'k-', 'LineWidth', 2);
set(h0, 'HandleVisibility','off');

set(gca, 'XTick', x0, 'XTickLabel', roi_list, 'XTickLabelRotation', 0);
xlabel('ROI');
ylabel('Mean slope (merged-voxel fit per subject)');
title('Mean slope values (MERGED voxels; Wilcoxon vs 0 & Con-Uncon, BH-FDR)');

% Axis appearance settings
ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'in';
ax.LineWidth = 1.5;
ax.FontSize  = 12;
grid off;

% Prevent clipping of the leftmost and rightmost bars
xlim([0.5, nRoi + 0.5]);

% y-axis limits (keep your original fixed range)
ylim([-1.2, 0.5]);

%% overlay subject dots
dot_color = [0 0 0];  % all subject dots are black
dot_size = 18;

if nSub == 1
    jitter = 0;
else
    jitter_span = 0.18;
    jitter = linspace(-jitter_span/2, jitter_span/2, nSub);
end

for s = 1:nSub
    % Con dots
    y = slope_con(s, :);
    m = isfinite(y);
    if any(m)
        scatter(xtips_con(m) + jitter(s), y(m), dot_size, 'o', ...
            'MarkerFaceColor', dot_color, 'MarkerEdgeColor', dot_color, 'LineWidth', 0.5, ...
            'HandleVisibility','off');
    end

    % Uncon dots
    y = slope_uncon(s, :);
    m = isfinite(y);
    if any(m)
        scatter(xtips_uncon(m) + jitter(s), y(m), dot_size, 'o', ...
            'MarkerFaceColor', dot_color, 'MarkerEdgeColor', dot_color, 'LineWidth', 0.5, ...
            'HandleVisibility','off');
    end
end

%% legend
legend([bh(1), bh(2)], {'Connected (Con)','Unconnected (Uncon)'}, ...
    'Location', 'bestoutside', 'Interpreter','none');

%% stars + lines: Con vs0 / Uncon vs0 / Con-Uncon
fprintf('\n=== FDR-corrected q values and star labels (fixed y=0.1; with lines) ===\n');

% ------- significance stars + lines (line width matches bar width) -------
y_star_fixed = 0.1;  % your original fixed y position for stars; here used as the line height
yL = ylim;
yR = yL(2) - yL(1);

% Place stars above the horizontal lines
y_sig_line = y_star_fixed;
y_sig_text = y_sig_line + 0.03 * yR;  % distance between stars and line
if y_sig_text > yL(2) - 0.01*yR
    y_sig_text = yL(2) - 0.01*yR;
end

halfLineW = barW/2;  % half-width of the horizontal line (same as bar width)

% Put the difference test (Con-Uncon) on a higher layer to avoid overlap with single-condition stars
y_sig_line_diff = y_sig_line + 0.06 * yR;
y_sig_text_diff = y_sig_line_diff + 0.03 * yR;
if y_sig_text_diff > yL(2) - 0.01*yR
    y_sig_text_diff = yL(2) - 0.01*yR;
end

for r = 1:nRoi
    % Con vs 0
    if isfinite(p_con_fdr(r))
        stars_con = p2stars(p_con_fdr(r)); % use q values for star labeling
        fprintf('ROI %s, Con vs0:     q = %.4g, stars = "%s"\n', roi_list{r}, p_con_fdr(r), stars_con);
        if ~isempty(stars_con)
            hL = line([xtips_con(r)-halfLineW, xtips_con(r)+halfLineW], [y_sig_line, y_sig_line], ...
                 'Color',[0 0 0], 'LineWidth', 1.5, 'Clipping','on');
            set(hL, 'HandleVisibility','off');

            text(xtips_con(r), y_sig_text, stars_con, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'Color',[0 0 0], 'FontSize', 18, 'FontWeight','bold', 'Clipping','on');
        end
    end

    % Uncon vs 0
    if isfinite(p_uncon_fdr(r))
        stars_uncon = p2stars(p_uncon_fdr(r));
        fprintf('ROI %s, Uncon vs0:   q = %.4g, stars = "%s"\n', roi_list{r}, p_uncon_fdr(r), stars_uncon);
        if ~isempty(stars_uncon)
            hL = line([xtips_uncon(r)-halfLineW, xtips_uncon(r)+halfLineW], [y_sig_line, y_sig_line], ...
                 'Color',[0 0 0], 'LineWidth', 1.5, 'Clipping','on');
            set(hL, 'HandleVisibility','off');

            text(xtips_uncon(r), y_sig_text, stars_uncon, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'Color',[0 0 0], 'FontSize', 18, 'FontWeight','bold', 'Clipping','on');
        end
    end

    % Con vs Uncon 
    if isfinite(p_diff_fdr(r))
        stars_diff = p2stars(p_diff_fdr(r));
        fprintf('ROI %s, Con-Uncon:   q = %.4g, stars = "%s"\n', roi_list{r}, p_diff_fdr(r), stars_diff);
        if ~isempty(stars_diff)
            x_mid = mean([xtips_con(r), xtips_uncon(r)]);

            hL = line([xtips_con(r), xtips_uncon(r)], [y_sig_line_diff, y_sig_line_diff], ...
                 'Color',[0 0 0], 'LineWidth', 1.5, 'Clipping','on');
            set(hL, 'HandleVisibility','off');

            text(x_mid, y_sig_text_diff, stars_diff, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'Color',[0 0 0], 'FontSize', 16, 'FontWeight','bold', 'Clipping','on');
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
    if nargin < 2, alpha = 0.05; end

    p = pvals(:);
    m = numel(p);
    [ps, idx] = sort(p);
    r = (1:m)';

    thresh = r/m * alpha;
    below  = find(ps <= thresh);

    h = false(m,1);
    if ~isempty(below)
        k = below(end);
        crit_p = ps(k);
        h(1:k) = true;
    else
        crit_p = 0;
    end

    adj_ps = ps .* m ./ r;
    adj_ps = flipud(cummin(flipud(adj_ps)));
    adj_ps(adj_ps>1) = 1;

    adj_p = nan(m,1);
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