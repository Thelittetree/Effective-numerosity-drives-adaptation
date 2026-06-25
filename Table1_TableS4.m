%% ============================================================
% Fit LME model only for Log DV from: T_diff_LME_withLinearAndLogDV.mat
%
% Model (voxel-level):
%   y = beta0 + beta_cond*x + map + x:map
%       + (1 + x | subject) + (1 | VoxelKey) + error
%
% Here:
%   y = y_log_PNdiff = log(PNadapt) - log(PNctrl)
%
% Notes:
% - x is coded as -0.5/+0.5
% - For fixed-effect model comparison, use FitMethod='ML'
%% ============================================================

% clear; clc;
matFile = 'T_diff_LME_withLinearAndLogDV.mat';
S = load(matFile);
assert(isfield(S,'T_lme') && istable(S.T_lme), 'do not found T_lme(table)。');
T = S.T_lme;

%% --------- Basic type hygiene ----------
if ~iscategorical(T.subject),   T.subject   = categorical(T.subject);   end
if ~iscategorical(T.roi),       T.roi       = categorical(T.roi);       end
if ~iscategorical(T.condition), T.condition = categorical(T.condition); end
if ~iscategorical(T.VoxelKey),  T.VoxelKey  = categorical(T.VoxelKey);  end

% Keep only con/uncon rows
T = T(ismember(string(T.condition), ["con","uncon"]), :);
T.condition = categorical(string(T.condition), ["con","uncon"]);

% Drop undefined ROI rows if any
mask_roi_ok = ~isundefined(T.roi);
if any(~mask_roi_ok)
    warning('Dropping %d rows with undefined ROI category.', sum(~mask_roi_ok));
    T = T(mask_roi_ok,:);
end

% Remove unused categories
T.subject   = removecats(T.subject);
T.roi       = removecats(T.roi);
T.condition = removecats(T.condition);
T.VoxelKey  = removecats(T.VoxelKey);

%% --------- Create centered condition regressor ----------
% con = +0.5, uncon = -0.5
x = nan(height(T),1);
x(T.condition=="con")   = +0.5;
x(T.condition=="uncon") = -0.5;
T.cond_c = x;

%% --------- Prepare Log DV ----------
Tlog = T;
mask_log = isfinite(Tlog.y_log_PNdiff);
Tlog = Tlog(mask_log,:);
Tlog.y = Tlog.y_log_PNdiff;

%% --------- Random-effects structure ----------
re_subject = '(1 + cond_c|subject)';
% If convergence/singularity occurs, use:
% re_subject = '(1|subject) + (cond_c-1|subject)';

re_voxel = '(1|VoxelKey)';

%% --------- Fixed-effects formulas ----------
f_full  = sprintf('y ~ cond_c*roi + %s + %s', re_subject, re_voxel);
f_noint = sprintf('y ~ cond_c + roi + %s + %s', re_subject, re_voxel);

%% ============================================================
% Fit model for Log DV
%% ============================================================
fprintf('\n===== LME: Log DV (logPNadapt - logPNctrl) =====\n');

lme_full_ML  = fitlme(Tlog, f_full,  'FitMethod','ML', 'DummyVarCoding','effects');
lme_noint_ML = fitlme(Tlog, f_noint, 'FitMethod','ML', 'DummyVarCoding','effects');

%% --------- ANOVA table for full model ----------
a = anova(lme_full_ML, 'DFMethod','satterthwaite');
disp('--- ANOVA (full model) ---');
disp(a);

%% --------- Likelihood-ratio test for interaction ----------
cmp = compare(lme_noint_ML, lme_full_ML);
disp('--- Likelihood-ratio test: add (cond_c:roi) ---');
disp(cmp);

%% --------- Optional: save results ----------
out_models = 'LME_fit_log_only.mat';
save(out_models, 'lme_full_ML', 'lme_noint_ML', 'cmp', 'a', '-v7.3');
fprintf('Saved model objects: %s\n', out_models);