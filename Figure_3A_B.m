%% =========================================================
% Scatter plots of preferred numerosity (connect / unconnect / control)
% across all ROIs, with percentage summaries for points above the identity line.
%
% This script:
% 1) Filters voxels using ve_control >= 0.3 and preferred numerosity values
%    within the valid range [1.05, 7.05] across all three conditions.
% 2) For each ROI, generates three scatter plots:
%       Fig1: Unconnect vs Control   (Y = unconnect, X = control, blue)
%       Fig2: Connect   vs Control   (Y = connect,   X = control, red)
%       Fig3: Unconnect vs Connect   (Y = unconnect, X = connect, gray)
% 3) Computes, for each ROI, the percentage of voxels with y > x in each
%    comparison and stores these values for later summary.
%% =========================================================
clear; clc;

data_dir = 'params\';
in_file  = fullfile(data_dir, 'sub02_merged_data.mat');
load(in_file, 'pref', 've', 'width'); 

roi_order = {'NTO','NPO','NPC1','NPC2','NPC3','NF'};
nROI = numel(roi_order);

prop_U_ctrl = nan(1, nROI);
prop_C_ctrl = nan(1, nROI);
prop_U_C    = nan(1, nROI);

xrange = [1 7];
yrange = [1 7];

for r = 1:nROI
    roi_name = roi_order{r};
    if ~isfield(pref, roi_name)
        warning('ROI %s not found in pref struct. Skip.', roi_name);
        continue;
    end

    fprintf('ROI %s\n', roi_name);

    pn_ctrl  = pref.(roi_name).control(:);
    pn_con   = pref.(roi_name).con(:);
    pn_uncon = pref.(roi_name).uncon(:);
    ve_ctrl  = ve.(roi_name).control(:);

    mask = ve_ctrl >= 0.3 & ...
           pn_ctrl  >= 1.05 & pn_ctrl  <= 7.05 & ...
           pn_con   >= 1.05 & pn_con   <= 7.05 & ...
           pn_uncon >= 1.05 & pn_uncon <= 7.05 & ...
           ~isnan(pn_ctrl)  & ~isnan(pn_con) & ~isnan(pn_uncon);

    pn_ctrl_f  = pn_ctrl(mask);
    pn_con_f   = pn_con(mask);
    pn_uncon_f = pn_uncon(mask);

    n_points = numel(pn_ctrl_f);
    if n_points == 0
        warning('ROI %s: no voxel passed the filter.', roi_name);
        continue;
    end

    figure;
    scatter(pn_ctrl_f, pn_uncon_f, 5, [0.07 0.62 1.00], 'filled');
    hold on;
    plot(xrange, yrange, 'k-');
    xlim(xrange); ylim(yrange);
    axis square;
    xlabel('Preferred Numerosity (Control condition)');
    ylabel('Preferred Numerosity (Unconnect condition)');
    title(sprintf('%s: Unconnect vs Control', roi_name));
    set(gca, 'Box','off', ...
             'TickDir','in', ...
             'LineWidth',1.5, ...
             'FontSize',12);

    pct_uncon_ctrl = mean(pn_uncon_f > pn_ctrl_f) * 100;
    fprintf('  Unconnect vs Control: %.2f %% points with y > x\n', pct_uncon_ctrl);
    prop_U_ctrl(r) = pct_uncon_ctrl;

    figure;
    scatter(pn_ctrl_f, pn_con_f, 5, [1 0 0], 'filled');
    hold on;
    plot(xrange, yrange, 'k-');
    xlim(xrange); ylim(yrange);
    axis square;
    xlabel('Preferred Numerosity (Control condition)');
    ylabel('Preferred Numerosity (Connect condition)');
    title(sprintf('%s: Connect vs Control', roi_name));
    set(gca, 'Box','off', ...
             'TickDir','in', ...
             'LineWidth',1.5, ...
             'FontSize',12);

    pct_con_ctrl = mean(pn_con_f > pn_ctrl_f) * 100;
    fprintf('  Connect   vs Control: %.2f %% points with y > x\n', pct_con_ctrl);
    prop_C_ctrl(r) = pct_con_ctrl;

    figure;
    scatter(pn_con_f, pn_uncon_f, 5, [0.5 0.5 0.5], 'filled');
    hold on;
    plot(xrange, yrange, 'k-');
    xlim(xrange); ylim(yrange);
    axis square;
    xlabel('Preferred Numerosity (Connect condition)');
    ylabel('Preferred Numerosity (Unconnect condition)');
    title(sprintf('%s: Unconnect vs Connect', roi_name));
    set(gca, 'Box','off', ...
             'TickDir','in', ...
             'LineWidth',1.5, ...
             'FontSize',12);

    pct_uncon_con = mean(pn_uncon_f > pn_con_f) * 100;
    fprintf('  Unconnect vs Connect: %.2f %% points with y > x\n\n', pct_uncon_con);
    prop_U_C(r) = pct_uncon_con;
end