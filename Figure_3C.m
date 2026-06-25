% Required existing variables:
% vals_delta_A, roi_cat_A, cond_cat_A, isOut_A, roi_order, nROI

cond_names_adapt = {'Unconnected','Connected'};
cond_colors = [0.07,0.62,1.00;
               1.00,0.00,0.00];
nCond_adapt = numel(cond_names_adapt);

figure('Name','Delta median PN','Color','w');

boxplot(vals_delta_A, {roi_cat_A, cond_cat_A}, ...
        'labelverbosity','all', ...
        'Symbol','o', ...
        'OutlierSize', 4);

xlabel('ROI', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Change in Preferred Numerosity', 'FontSize', 14, 'FontWeight', 'bold');

grid off;
set(gca, 'YGrid','off','XGrid','off');
hold on;

yline(0, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'HandleVisibility','off');

hBox = findobj(gca, 'Tag', 'Box');
nBox = numel(hBox);

boxCenters = nan(nROI, nCond_adapt);
h_uncon = [];
h_con   = [];

for r = 1:nROI
    for c = 1:nCond_adapt
        idx = (r-1)*nCond_adapt + c;
        if idx > nBox, continue; end
        hb = hBox(nBox - idx + 1);
        thisX = get(hb, 'XData');
        thisY = get(hb, 'YData');

        boxCenters(r,c) = mean(thisX(:));

        thisColor = cond_colors(c,:);
        hp = patch(thisX, thisY, thisColor, ...
            'FaceAlpha', 0.4, ...
            'EdgeColor', thisColor);

        if r == 1
            if c == 1 && isempty(h_uncon), h_uncon = hp; end
            if c == 2 && isempty(h_con),   h_con   = hp; end
        end
    end
end

ROI_centers = mean(boxCenters, 2, 'omitnan');
set(gca, 'XTick', ROI_centers, 'XTickLabel', roi_order);

set(gca, 'Box','off', 'TickDir','in', 'LineWidth',1.5, 'FontSize',12);

if ~isempty(h_uncon) && ~isempty(h_con)
    legend([h_uncon, h_con], ...
        {'Unconnected Adaptation Condition', 'Connected Adaptation Condition'}, ...
        'Location','best', 'Box','off');
end
hold off;