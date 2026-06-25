%% Data analysis of behavioral experiment
% Version: integrated and corrected
% Full dataset: 13 × 4
% Columns:
%   1 = No adaptation
%   2 = Connected adaptation
%   3 = Unconnected adaptation
%   4 = PSE in the "40 connected dots" condition
%
% Functions:
% 1) One-sided Wilcoxon signed-rank test against fixed reference value 40
%    using underestimation = 40 - PSE
%    (alternative hypothesis: 40 - PSE > 0, i.e., PSE < 40)
%    + Hodges-Lehmann estimate of underestimation (40 - PSE)
%    + 95% CI
% 2) Shapiro-Wilk test (built-in Monte Carlo version, optional)
% 3) Mauchly's test of sphericity
% 4) rm-ANOVA main effect of condition
% 5) Paired t-tests + Holm correction
% 6) Plot Figure 1/2

clear; clc; close all;

%% ================== Parameters (editable) ==================
ALPHA_WILCOXON   = 0.05;   % one-sided Wilcoxon alpha
ALPHA_CI         = 0.05;   % 95% CI for HL estimate

ALPHA_NORMALITY  = 0.05;
ALPHA_SPHERICITY = 0.05;
ALPHA_AOV        = 0.05;

DO_SHAPIRO_WILK = true;  % true: run built-in Monte Carlo SW (slower); false: skip
SW_SEED = 1;
SW_WSIM = 30000;         
SW_PSIM = 20000;         

REF_NUM = 40;            % fixed physical numerosity reference

%% ================== Data ==================
PSE_data_all = [ ...
    17.99 18.80 13.36 30;
    17.13 15.52 17.17 30;
    17.78 15.66 14.30 29;
    18.32 15.93 16.00 30;
    19.71 15.89 14.86 27;
    19.20 20.00 17.16 20;
    15.17 16.73 14.57 20;
    16.12 13.63 14.75 30;
    17.06 16.95 14.92 30;
    18.43 16.27 14.15 30;
    18.26 16.81 12.24 30;
    16.00 15.18 15.84 19;
    16.15 13.98 11.28 30
];

% First 3 columns are used for rm-ANOVA / post hoc / plots
PSE_data = PSE_data_all(:,1:3);

% Last column is used for the one-sample Wilcoxon signed-rank test against 40
PSE_40connected = PSE_data_all(:,4);

conditions_short = {'NoAdapt','Connected','Unconnected'};
cond_labels_en   = {'No adaptation','Connected adaptation','Unconnected adaptation'};

nSub = size(PSE_data_all,1);
subjects = (1:nSub)';

%% 1) One-sided Wilcoxon signed-rank test against fixed reference 40
% Underestimation: positive values indicate PSE < 40
underestimation40 = REF_NUM - PSE_40connected;

% One-sided Wilcoxon signed-rank test:
% H1: underestimation > 0  <=>  40 - PSE > 0  <=>  PSE < 40
[p_wil, h_wil, stats_wil] = signrank(underestimation40, 0, 'tail', 'right', 'alpha', ALPHA_WILCOXON);

if isfield(stats_wil, 'signedrank')
    V_wil = stats_wil.signedrank;
else
    V_wil = NaN;
end

% Descriptives for the PSE values in the 40-connected-dots condition
mean_pse40 = mean(PSE_40connected, 'omitnan');
sd_pse40   = std(PSE_40connected, 0, 'omitnan');   % sample SD

% Descriptives for underestimation
median_under_raw = median(underestimation40, 'omitnan');

% Hodges-Lehmann estimate and its 95% CI for underestimation = 40 - PSE
[hl_est_under, hl_ci_under] = hodges_lehmann_1sample_ci(underestimation40, ALPHA_CI);

%% 2) rm-ANOVA preparation (first 3 columns only)
T = array2table(PSE_data, 'VariableNames', conditions_short);
T.Subject = subjects;

within = table(categorical(conditions_short(:), conditions_short, 'Ordinal', true), ...
    'VariableNames', {'Condition'});

rm = fitrm(T, 'NoAdapt-Unconnected ~ 1', 'WithinDesign', within);
aovTbl = ranova(rm, 'WithinModel', 'Condition');

%% Pairwise comparisons (paired t-tests + Holm correction)
pairs_list = { ...
    'NoAdapt','Connected'; ...
    'NoAdapt','Unconnected'; ...
    'Connected','Unconnected' };

p_unc = nan(3,1);
for k = 1:3
    A = pairs_list{k,1};
    B = pairs_list{k,2};
    [~, p] = ttest(T.(A), T.(B)); % paired t-test
    p_unc(k) = p;
end
p_corr = holm_adjust(p_unc);

% Corresponding to: 1-2, 1-3, 2-3
p_12 = p_corr(1);
p_13 = p_corr(2);
p_23 = p_corr(3);

%% Descriptive statistics (first 3 columns)
means = mean(PSE_data, 1, 'omitnan');
sds   = std(PSE_data,  0, 1, 'omitnan');

%% 3) Normality (Shapiro-Wilk) + sphericity (Mauchly)
sw_W = nan(1,3);
sw_p = nan(1,3);
if DO_SHAPIRO_WILK
    for j = 1:3
        xj = PSE_data(:,j);
        [sw_W(j), sw_p(j)] = shapiro_wilk_mc(xj, SW_WSIM, SW_PSIM, SW_SEED + j);
    end
end
normal_ok = all(sw_p > ALPHA_NORMALITY);

try
    mauchlyTbl = mauchly(rm);
    chi2_m = mauchlyTbl.ChiStat(1);
    df_m   = mauchlyTbl.DF(1);
    p_m    = mauchlyTbl.pValue(1);
catch
    chi2_m = NaN; df_m = NaN; p_m = NaN;
end
sphericity_ok = (~isnan(p_m)) && (p_m > ALPHA_SPHERICITY);


%% 4) rm-ANOVA main effect
rn = aovTbl.Properties.RowNames;

% Main effect row: usually "(Intercept):Condition"
iEff = find(contains(rn,'Condition') & ~contains(rn,'Error'), 1);

% Error row: usually "Error(Condition)"
iErr = find(contains(rn,'Error(Condition)'), 1);

if isempty(iEff) || isempty(iErr)
    error('Unable to automatically locate the Condition main-effect row or Error(Condition) row in the ranova table. Please check aovTbl row names.');
end

F_val = aovTbl.F(iEff);
df1   = aovTbl.DF(iEff);
df2   = aovTbl.DF(iErr);

% Sphericity handling: if violated -> use GG corrected p
if ~isnan(p_m) && p_m < ALPHA_SPHERICITY
    p_aov = aovTbl.pValueGG(iEff);
else
    p_aov = aovTbl.pValue(iEff);
end

SS_eff = aovTbl.SumSq(iEff);
SS_err = aovTbl.SumSq(iErr);
eta_p2 = SS_eff / (SS_eff + SS_err);

main_effect_sig = (p_aov < ALPHA_AOV);

%% Printed output
fprintf('\n==================== SUMMARY (for writing) ====================\n');
fprintf('N subjects = %d\n', nSub);

% ---------- Fixed-reference Wilcoxon test ----------
fprintf('\n[Single-sample test against %g] One-sided Wilcoxon signed-rank test\n', REF_NUM);
fprintf('  Underestimation variable: %g - PSE\n', REF_NUM);
fprintf('  Alternative hypothesis: %g - PSE > 0 (equivalent to PSE < %g)\n', REF_NUM, REF_NUM);
fprintf('  PSEs in the 40-connected-dots condition: mean = %.2f, SD = %.2f\n', mean_pse40, sd_pse40);
fprintf('  Raw median underestimation (%g - PSE) = %.2f\n', REF_NUM, median_under_raw);
fprintf('  Hodges-Lehmann estimate of underestimation = %.2f\n', hl_est_under);
fprintf('  95%% CI of HL estimate = [%.2f, %.2f]\n', hl_ci_under(1), hl_ci_under(2));
fprintf('  One-sided Wilcoxon signed-rank test: V = %.0f, p %s\n', V_wil, p_text_with_symbol(p_wil));
fprintf('  ==> Significant underestimation? %s\n', yesno(h_wil == 1));

% ---------- Auto-generated English sentence for the first part ----------
fprintf('\n[English sentence for Results: fixed-reference test]\n');
fprintf(['Across participants, the unconnected array that matched the perceived numerosity of %g connected dots ' ...
         'yielded PSEs below %g (mean = %.2f, SD = %.2f), corresponding to a median underestimation of %.2f dots ' ...
         '(Hodges-Lehmann estimate = %.2f, 95%% CI [%.2f, %.2f]). This underestimation was significant ' ...
         '(one-sided Wilcoxon signed-rank test: V = %.0f, p %s), establishing a dissociation between physical and perceived numerosity induced by connectedness.\n'], ...
         REF_NUM, REF_NUM, mean_pse40, sd_pse40, median_under_raw, hl_est_under, hl_ci_under(1), hl_ci_under(2), V_wil, p_text_with_symbol(p_wil));

% ---------- Normality ----------
if DO_SHAPIRO_WILK
    fprintf('\n[Normality] Shapiro-Wilk (per condition)\n');
    for j = 1:3
        fprintf('  %s: W = %.4f, p = %.4g\n', conditions_short{j}, sw_W(j), sw_p(j));
    end
    fprintf('  ==> All p > %.3f ? %s\n', ALPHA_NORMALITY, yesno(normal_ok));
else
    fprintf('\n[Normality] Shapiro-Wilk was skipped (DO_SHAPIRO_WILK = false).\n');
end

% ---------- Sphericity ----------
fprintf('\n[Sphericity] Mauchly''s test\n');
if ~isnan(p_m)
    fprintf('  Chi-square(%d) = %.2f, p = %.4g  ==> Not violated (p > %.3f)? %s\n', ...
        df_m, chi2_m, p_m, ALPHA_SPHERICITY, yesno(sphericity_ok));
else
    fprintf('  mauchly(rm) not available. (Chi-square / df / p not computed)\n');
end

% ---------- rm-ANOVA ----------
fprintf('\n[rm-ANOVA] Main effect of Condition\n');
fprintf('  F(%d, %d) = %.3f, p %s, partial eta^2 = %.2f\n', df1, df2, F_val, p_text_with_symbol(p_aov), eta_p2);

% ---------- Descriptives ----------
fprintf('\n[Descriptives] Mean ± SD\n');
fprintf('  NoAdapt     : Mean = %.2f, SD = %.2f\n', means(1), sds(1));
fprintf('  Connected   : Mean = %.2f, SD = %.2f\n', means(2), sds(2));
fprintf('  Unconnected : Mean = %.2f, SD = %.2f\n', means(3), sds(3));

% ---------- Pairwise ----------
fprintf('\n[Pairwise] Paired t-tests with Holm correction\n');
fprintf('  NoAdapt vs Connected     : p %s\n', p_text_with_symbol(p_12));
fprintf('  NoAdapt vs Unconnected   : p %s\n', p_text_with_symbol(p_13));
fprintf('  Connected vs Unconnected : p %s\n', p_text_with_symbol(p_23));

%% ============================================================
%% Figure 1
%% ============================================================
ref_y = 15.0;
fig1 = figure('Color','w'); ax1 = axes(fig1); hold(ax1,'on');

palette = lines(10);
for s = 1:nSub
    y = PSE_data(s,:);
    c = palette(mod(s-1, size(palette,1)) + 1, :);
    plot(ax1, [1 2 3], y, ':o', 'Color', c, 'MarkerSize', 4, 'LineWidth', 1.0);
end

yline(ax1, ref_y, '--', 'LineWidth', 1.2);

xlim(ax1, [0.6 3.4]);
ylim(ax1, [10 max(22, max(PSE_data(:)) + 2.5)]);

xticks(ax1, [1 2 3]);
xticklabels(ax1, cond_labels_en);

xlabel(ax1, 'Condition', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax1, 'PSE value (dots)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, 'PSE across adaptation conditions (individual trajectories)', 'FontSize', 13, 'FontWeight', 'bold');

set(ax1, 'Box','off', 'TickDir','in', 'LineWidth',1.5, 'FontSize',12);

%% ============================================================
%% Figure 2
%% ============================================================
fig2 = figure('Color','w'); ax2 = axes(fig2); hold(ax2,'on');

x = 1:3;

barColor = [0.30 0.60 0.90];
dotColor = [0.20 0.20 0.20];
rng(0);
jitter = 0.06;

sig_y0  = 21;
sig_h   = 0.18;
sig_gap = 0.55;

xlim(ax2, [0.4 3.6]);
ylim(ax2, [12 24]);

hb = bar(ax2, x, means, 0.65, ...
    'FaceColor', barColor, 'EdgeColor','k', 'LineWidth', 0.8);
hb.DisplayName = 'Mean';

he = errorbar(ax2, x, means, sds, ...
    'LineStyle', 'none', 'Color', 'r', 'LineWidth', 2, 'CapSize', 10);
he.DisplayName = 'Standard deviation';

hs_first = [];
for j = 1:3
    yj = PSE_data(:,j);
    xj = x(j) + (rand(nSub,1)-0.5)*2*jitter;
    hs = scatter(ax2, xj, yj, 18, dotColor, 'filled');
    if isempty(hs_first)
        hs_first = hs;
        hs_first.DisplayName = 'Individual data points';
    else
        set(hs, 'HandleVisibility','off');
    end
end

yline(ax2, ref_y, '--', 'LineWidth', 1.2);

xticks(ax2, x);
xticklabels(ax2, cond_labels_en);

xlabel(ax2, 'Condition', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax2, 'Mean PSE value (dots)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax2, 'PSE across adaptation conditions', 'FontSize', 13, 'FontWeight', 'bold');

lgd = legend(ax2, [hb he hs_first], {'Mean','Standard deviation','Individual data points'}, ...
    'Orientation','horizontal');
set(lgd, 'Box','on');
lgd.AutoUpdate = 'off';

lgd.Units = 'normalized';
axPos  = ax2.Position;
lgdPos = lgd.Position;
lgdPos(1) = axPos(1) + (axPos(3) - lgdPos(3)) / 2;
lgdPos(2) = axPos(2) + axPos(4) - lgdPos(4) - 0.005;
lgd.Position = lgdPos;

set(ax2, 'Box','off', 'TickDir','in', 'LineWidth',1.5, 'FontSize',12);

% 1-2 and 2-3 are drawn at the same height
y12 = sig_y0;
y23 = sig_y0;
y13 = sig_y0 + sig_gap;

add_sig_bracket(ax2, 1, 2, y12, sig_h, p_12);
add_sig_bracket(ax2, 2, 3, y23, sig_h, p_23);
add_sig_bracket(ax2, 1, 3, y13, sig_h, p_13);

%% ======================= Local functions =======================
function s = yesno(tf)
    if tf, s = 'YES'; else, s = 'NO'; end
end

function s = p_text_with_symbol(p)
% Return "< 0.001" or "= 0.047" (used as "p %s")
    if isnan(p)
        s = 'NA';
        return;
    end
    if p < 0.001
        s = '< 0.001';
    else
        s = sprintf('= %.3f', p);
    end
end

function s = p_to_star_en(p)
    if isnan(p), s = 'n.s.'; return; end
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

function add_sig_bracket(ax, x1, x2, y, h, p)
    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k', ...
        'LineWidth',1.0, 'Clipping','off', 'HandleVisibility','off');
    text(ax, (x1+x2)/2, y+h, p_to_star_en(p), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',12, 'FontWeight','bold', 'HandleVisibility','off');
end

function p_adj = holm_adjust(p)
    p = p(:);
    m = numel(p);
    [ps, idx] = sort(p, 'ascend');

    adj = nan(m,1);
    for i = 1:m
        adj(i) = (m - i + 1) * ps(i);
    end
    for i = 2:m
        adj(i) = max(adj(i), adj(i-1));
    end
    adj = min(adj, 1);

    p_adj = nan(m,1);
    p_adj(idx) = adj;
end

function [hl_est, ci] = hodges_lehmann_1sample_ci(d, alpha)
% Hodges-Lehmann estimate and 100*(1-alpha)% CI for one-sample data d
% Here d = underestimation = 40 - PSE
%
% HL estimate:
%   median of Walsh averages: (d_i + d_j)/2, i <= j
%
% CI:
%   obtained by inverting the two-sided Wilcoxon signed-rank test
%   over the sorted unique Walsh averages.

    d = d(:);
    d = d(~isnan(d));

    n = numel(d);
    if n < 1
        hl_est = NaN;
        ci = [NaN NaN];
        return;
    end

    % Walsh averages
    Nw = n*(n+1)/2;
    walsh = nan(Nw,1);
    idx = 1;
    for i = 1:n
        for j = i:n
            walsh(idx) = (d(i) + d(j)) / 2;
            idx = idx + 1;
        end
    end
    walsh = sort(walsh);

    % HL estimate
    hl_est = median(walsh);

    % CI by inversion of the two-sided signrank test
    cand = unique(walsh);
    pvals = nan(size(cand));
    for k = 1:numel(cand)
        pvals(k) = signrank(d, cand(k), 'tail', 'both');
    end

    keep = (pvals >= alpha);
    if any(keep)
        ci = [cand(find(keep,1,'first')), cand(find(keep,1,'last'))];
    else
        ci = [NaN NaN];
    end
end

function [Wobs, pval] = shapiro_wilk_mc(x, nSimW, nSimP, seed)
% Monte Carlo version of Shapiro-Wilk (self-contained)
% Left-tailed: smaller W indicates stronger deviation from normality
% p = P(Wsim <= Wobs)

    x = x(:);
    x = x(~isnan(x));
    n = numel(x);
    if n < 3
        Wobs = NaN; pval = NaN; return;
    end

    rng(seed);

    % (1) Estimate m and V
    Z = randn(nSimW, n);
    Zs = sort(Z, 2);
    m  = mean(Zs, 1)';      % n×1
    V  = cov(Zs);           % n×n

    % (2) Weights a
    b = V \ m;
    a = b / norm(b);

    % (3) Observed W
    xs = sort(x);
    xbar = mean(xs);
    denom = sum((xs - xbar).^2);
    numer = (a' * xs)^2;
    Wobs = numer / denom;

    % (4) Simulate null distribution to estimate p
    Z2 = randn(nSimP, n);
    Z2s = sort(Z2, 2);
    Z2bar = mean(Z2s, 2);
    denom2 = sum((Z2s - Z2bar).^2, 2);
    numer2 = (Z2s * a).^2;
    Wsim = numer2 ./ denom2;

    pval = (sum(Wsim <= Wobs) + 1) / (nSimP + 1);
end