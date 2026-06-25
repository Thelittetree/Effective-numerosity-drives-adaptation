%% Figure S1: Additional behavioral experiment and results
% rows = subjects, columns = conditions
% col order:
% 1 no adapter | 2 contour | 3 connected | 4 middle | 5 unconnected | 6 lines
PSE_data = [
    18.2300   15.8100   19.2100   15.0200   17.8400   16.8300;
    20.0000   20.0000   17.2500   17.4500   17.1300   17.5700;
    16.8200   17.5300   17.2800   15.7300   14.8600   16.2400;
    17.4500   19.7200   20.0000   20.0000   20.0000   18.6000;
    16.5900   14.6100   16.8800   15.9400   14.8600   14.5300;
    17.8600   17.5900   16.9100   15.1500   16.5000   16.1600;
    15.5600   16.5500   14.5700   14.8800   16.6800   15.2300;
    17.5200   14.5500   15.4600   16.3000   16.1200   16.1700;
    19.1400   16.1100   16.1500   19.8500   18.1800   18.1300;
    17.1200   19.2000   20.0000   19.7500   17.2100   16.1100;
    20.0000   18.9500   18.2600   15.6300   19.3400   15.5700;
    14.4200   16.5700   16.0700   14.8200   17.6400   13.3400;
    18.3900   16.2500   13.9900   17.4300   15.8500   15.4400;
    16.1500   13.1000   13.9800   14.4000   11.2800   11.1200;
    17.9900   17.4000   18.8000   14.4900   13.3600   16.3600;
    17.1300   15.7500   15.5200   14.1400   17.1700   18.6600;
    17.7800   16.9200   15.6600   15.1600   14.3000   14.4000;
    18.3200   17.9500   15.9300   16.8700   16.0000   18.8400;
    19.7100   17.7200   15.8900   14.8800   14.8600   18.9000;
    19.2000   19.3500   20.0000   18.7700   17.1600   18.8300;
    16.1200   13.7100   13.6300   17.8400   14.7500   16.0200;
    17.0600   17.5600   16.9500   17.5100   14.9200   19.2800;
    18.4300   12.2300   16.2700   14.1800   14.1500   12.1300;
    18.2600   11.9300   16.8100   10.6400   12.2400   11.8200;
    16.0000   18.1200   15.1800   16.9800   15.8400   17.6100;
    17.6400   19.2800   19.6800   19.0500   18.5200   18.0700;
    16.3800   14.9600   14.9500   12.3300   14.2100   13.6000
];

alpha = 0.05;

%% Extract four conditions (order: contour, connected, lines, unconnected)
x_contour     = PSE_data(:, 2);
x_connected   = PSE_data(:, 3);
x_lines       = PSE_data(:, 6);
x_unconnected = PSE_data(:, 5);

X = [x_contour, x_connected, x_lines, x_unconnected];

% Merge groups (mean within each subject)
grp_CC = mean([x_connected, x_contour], 2);     % connected + contour
grp_LU = mean([x_lines, x_unconnected], 2);     % lines + unconnected

%% Color settings
col_uncon   = [0.3, 0.6, 0.9];
col_con     = [0.3, 0.6, 0.9];
col_line    = [0.3, 0.6, 0.9];
col_contour = [0.3, 0.6, 0.9];

bar_colors = [
    col_contour;
    col_con;
    col_line;
    col_uncon
];


xpos = [1, 2.4, 6.0, 7.4];
mean_PSE = mean(X, 1, 'omitnan');
n_eff    = sum(isfinite(X), 1);
sem_PSE  = std(X, 0, 1, 'omitnan') ./ sqrt(n_eff);

%% Plot
figure;
hold on;

% Bar plot (wider BarWidth: 0.90)
barW = 0.90;
b = bar(xpos, mean_PSE, 'BarWidth', barW, 'FaceColor', 'flat');
b.CData = bar_colors;

% SEM error bars: red points/lines (others unchanged)
errorbar(xpos, mean_PSE, sem_PSE, 'r', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 15);

% Subject-level dots (without jitter)
num_subjects   = size(X, 1);
num_conditions = size(X, 2);
for i = 1:num_conditions
    scatter(ones(num_subjects, 1) * xpos(i), X(:, i), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Statistics: normality test (difference scores) -> t-test or signrank
warnState = warning; warning('off','all');

% (1) contour vs connected
mask12 = isfinite(x_contour) & isfinite(x_connected);
a1 = x_contour(mask12);
b1v = x_connected(mask12);
d1 = b1v - a1;  % connected - contour
[hL1, pL1] = lillietest(d1, 'Alpha', alpha);
[hJ1, pJ1] = jbtest(d1, alpha);
isNormal1 = (hL1==0) && (hJ1==0);
if isNormal1
    [~, p1] = ttest(d1, 0, 'Alpha', alpha);
    testName1 = 'paired t-test';
else
    p1 = signrank(a1, b1v, 'Alpha', alpha);
    testName1 = 'signrank';
end

% (2) lines vs unconnected
mask34 = isfinite(x_lines) & isfinite(x_unconnected);
a2 = x_lines(mask34);
b2v = x_unconnected(mask34);
d2 = b2v - a2;  % unconnected - lines
[hL2, pL2] = lillietest(d2, 'Alpha', alpha);
[hJ2, pJ2] = jbtest(d2, alpha);
isNormal2 = (hL2==0) && (hJ2==0);
if isNormal2
    [~, p2] = ttest(d2, 0, 'Alpha', alpha);
    testName2 = 'paired t-test';
else
    p2 = signrank(a2, b2v, 'Alpha', alpha);
    testName2 = 'signrank';
end

% (3) merged groups: mean(CC) vs mean(LU)
maskG = isfinite(grp_CC) & isfinite(grp_LU);
g1 = grp_CC(maskG);
g2 = grp_LU(maskG);
dG = g2 - g1;  % LU - CC
[hLG, pLG] = lillietest(dG, 'Alpha', alpha);
[hJG, pJG] = jbtest(dG, alpha);
isNormalG = (hLG==0) && (hJG==0);
if isNormalG
    [~, pG] = ttest(dG, 0, 'Alpha', alpha);
    testNameG = 'paired t-test';
else
    pG = signrank(g1, g2, 'Alpha', alpha);
    testNameG = 'signrank';
end

warning(warnState);

% Significance labels (non-significant -> n.s.)
if p1 < 1e-4, lab1 = '****';
elseif p1 < 1e-3, lab1 = '***';
elseif p1 < 1e-2, lab1 = '**';
elseif p1 < 0.05, lab1 = '*';
else, lab1 = 'n.s.'; end

if p2 < 1e-4, lab2 = '****';
elseif p2 < 1e-3, lab2 = '***';
elseif p2 < 1e-2, lab2 = '**';
elseif p2 < 0.05, lab2 = '*';
else, lab2 = 'n.s.'; end

if pG < 1e-4, labG = '****';
elseif pG < 1e-3, labG = '***';
elseif pG < 1e-2, labG = '**';
elseif pG < 0.05, labG = '*';
else, labG = 'n.s.'; end

% Bracket drawing
ylim([10 21]);

% Draw dashed reference line at y = 15
yline(15, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'HandleVisibility','off');

dh  = 0.12;
gap = 0.1;
lw1 = 1.8;
lw2 = 2.2;

yMaxData = max(X(:), [], 'omitnan');
yMaxBar  = max(mean_PSE + sem_PSE);
baseTop  = max([yMaxData, yMaxBar]);

y12    = baseTop + 0.18;
yMerge = y12 + 0.30;
yTop   = yMerge + 0.30;

yNeededTop = yTop + dh + gap + 0.10;
if yNeededTop > 21
    shiftDown = yNeededTop - 21;
    y12    = y12    - shiftDown;
    yMerge = yMerge - shiftDown;
    yTop   = yTop   - shiftDown;
end

% (1) contour vs connected
plot([xpos(1) xpos(1) xpos(2) xpos(2)], [y12 y12+dh y12+dh y12], 'k', 'LineWidth', lw1);
text(mean([xpos(1) xpos(2)]), y12+dh+gap, lab1, 'HorizontalAlignment','center', 'FontSize', 12);

% (2) lines vs unconnected
plot([xpos(3) xpos(3) xpos(4) xpos(4)], [y12 y12+dh y12+dh y12], 'k', 'LineWidth', lw1);
text(mean([xpos(3) xpos(4)]), y12+dh+gap, lab2, 'HorizontalAlignment','center', 'FontSize', 12);

% (3) Summary brackets (without text)
plot([xpos(1) xpos(1) xpos(2) xpos(2)], [yMerge yMerge+dh yMerge+dh yMerge], 'k', 'LineWidth', lw2);
plot([xpos(3) xpos(3) xpos(4) xpos(4)], [yMerge yMerge+dh yMerge+dh yMerge], 'k', 'LineWidth', lw2);

% (4) Comparison across the two merged groups
plot([xpos(1) xpos(1) xpos(4) xpos(4)], [yTop yTop+dh yTop+dh yTop], 'k', 'LineWidth', lw2);
text(mean([xpos(1) xpos(4)]), yTop+dh+gap, labG, 'HorizontalAlignment','center', 'FontSize', 12);

% Aesthetic settings
set(gca, 'XTick', xpos, 'XTickLabel', {'contour','connected','lines','unconnected'});
xlabel('Adaptation Condition', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('PSE', 'FontSize', 12, 'FontWeight', 'bold');
title('PSE under Different Adaptation Conditions', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 11);
ax = gca;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;

grid off;

xlim([min(xpos)-0.8, max(xpos)+0.8]);
hold off;

%% Console output
fprintf('\n===== Stats summary (alpha=%.3f) =====\n', alpha);
fprintf('[contour vs connected]  test=%s | p=%.6g | n=%d\n', testName1, p1, numel(d1));
fprintf('[lines vs unconnected] test=%s | p=%.6g | n=%d\n', testName2, p2, numel(d2));
fprintf('[merged CC vs merged LU] test=%s | p=%.6g | n=%d\n', testNameG, pG, numel(dG));
fprintf('Normality p-values (Lillie/JB): pair1(%.3g/%.3g), pair2(%.3g/%.3g), merged(%.3g/%.3g)\n', ...
    pL1, pJ1, pL2, pJ2, pLG, pJG);