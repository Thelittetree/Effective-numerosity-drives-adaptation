%% ============================================================
% Fit voxel-level LME model for log-scale adaptation effect
%
% This script includes:
%   1. Loading the LME-ready table from:
%        data/params/T_diff_LME_withLinearAndLogDV.mat
%
%   2. Preparing the log-scale dependent variable:
%        y = y_log_PNdiff
%          = log(PNadapt) - log(PNctrl)
%
%   3. Fitting voxel-level linear mixed-effects models:
%
%        Full model:
%        y ~ cond_c * roi + (1 + cond_c | subject) + (1 | VoxelKey)
%
%        Reduced model without condition-by-map interaction:
%        y ~ cond_c + roi + (1 + cond_c | subject) + (1 | VoxelKey)
%
%   4. Testing fixed effects in the full model using ANOVA with
%      Satterthwaite degrees of freedom.
%
%   5. Comparing the full model against the reduced model using
%      likelihood-ratio test to evaluate the condition-by-map interaction.
%
%   6. Saving fitted model objects and statistical outputs to:
%        results/LME_fit_log_only.mat
%
% Notes:
%   - The analysis is performed at the voxel level.
%   - roi is used as the map factor.
%   - condition is restricted to con and uncon.
%   - cond_c is effect-coded:
%         con   = +0.5
%         uncon = -0.5
%   - For fixed-effect model comparison, models are fitted with
%     FitMethod = 'ML'.
%
% Expected project structure:
%
%   project_root/
%   |-- Fit_LME_log_only.m
%   |-- data/
%   |   |-- params/
%   |   |   |-- T_diff_LME_withLinearAndLogDV.mat
%   |-- results/
%
% If this script is placed inside an analysis/ subfolder, replace:
%     project_root = script_dir;
% with:
%     project_root = fileparts(script_dir);
%% ============================================================

clear; clc; close all;

% Paths
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

matFile = fullfile(data_dir, 'T_diff_LME_withLinearAndLogDV.mat');
out_models = fullfile(out_dir, 'LME_fit_log_only.mat');

if ~exist(matFile, 'file')
    error('Input file was not found:\n%s\nPlease check project_root and data folder structure.', matFile);
end


S = load(matFile);
assert(isfield(S, 'T_lme') && istable(S.T_lme), ...
    'T_lme table was not found in the input MAT file.');

T = S.T_lme;

fprintf('\nLoaded LME table from:\n%s\n', matFile);
fprintf('Original table size: %d rows × %d columns\n', height(T), width(T));

if ~iscategorical(T.subject)
    T.subject = categorical(T.subject);
end

if ~iscategorical(T.roi)
    T.roi = categorical(T.roi);
end

if ~iscategorical(T.condition)
    T.condition = categorical(T.condition);
end

if ~iscategorical(T.VoxelKey)
    T.VoxelKey = categorical(T.VoxelKey);
end

T = T(ismember(string(T.condition), ["con", "uncon"]), :);
T.condition = categorical(string(T.condition), ["con", "uncon"]);

fprintf('After keeping con/uncon rows: %d rows\n', height(T));
mask_roi_ok = ~isundefined(T.roi);

if any(~mask_roi_ok)
    warning('Dropping %d rows with undefined ROI category.', sum(~mask_roi_ok));
    T = T(mask_roi_ok, :);
end

T.subject   = removecats(T.subject);
T.roi       = removecats(T.roi);
T.condition = removecats(T.condition);
T.VoxelKey  = removecats(T.VoxelKey);

% Effect coding:
%   con   = +0.5
%   uncon = -0.5

x = nan(height(T), 1);

x(T.condition == "con")   = +0.5;
x(T.condition == "uncon") = -0.5;

T.cond_c = x;

if any(~isfinite(T.cond_c))
    error('Some rows have invalid cond_c values. Please check condition labels.');
end

%% Prepare log-scale dependent variable
% y_log_PNdiff = log(PNadapt) - log(PNctrl)

Tlog = T;

if ~ismember('y_log_PNdiff', Tlog.Properties.VariableNames)
    error('Variable y_log_PNdiff was not found in T_lme.');
end

mask_log = isfinite(Tlog.y_log_PNdiff);
n_removed_log = sum(~mask_log);

if n_removed_log > 0
    warning('Removing %d rows with non-finite y_log_PNdiff.', n_removed_log);
end

Tlog = Tlog(mask_log, :);
Tlog.y = Tlog.y_log_PNdiff;

fprintf('Final Log-DV table size: %d rows × %d columns\n', height(Tlog), width(Tlog));

% Random-effects structure
% Main random-effects structure:
%   subject-level random intercept and random condition slope
%   voxel-level random intercept

re_subject = '(1 + cond_c|subject)';

% If convergence or singularity occurs, use this alternative:
% re_subject = '(1|subject) + (cond_c-1|subject)';

re_voxel = '(1|VoxelKey)';

% Fixed-effects formulas-
% Full model:
%   condition × ROI interaction
%
% Reduced model:
%   condition main effect + ROI main effect only

f_full  = sprintf('y ~ cond_c*roi + %s + %s', re_subject, re_voxel);
f_noint = sprintf('y ~ cond_c + roi + %s + %s', re_subject, re_voxel);

fprintf('\nFull model formula:\n%s\n', f_full);
fprintf('\nReduced model formula:\n%s\n', f_noint);

% Fit LME models for Log DV

fprintf('\n===== LME: Log DV: log(PNadapt) - log(PNctrl) =====\n');

lme_full_ML = fitlme( ...
    Tlog, ...
    f_full, ...
    'FitMethod', 'ML', ...
    'DummyVarCoding', 'effects');

lme_noint_ML = fitlme( ...
    Tlog, ...
    f_noint, ...
    'FitMethod', 'ML', ...
    'DummyVarCoding', 'effects');

%% ANOVA table for the full model
fprintf('\n--- ANOVA for full model ---\n');

a = anova(lme_full_ML, 'DFMethod', 'satterthwaite');
disp(a);

%%  Likelihood-ratio test for interaction
fprintf('\n--- Likelihood-ratio test: add cond_c:roi interaction ---\n');

cmp = compare(lme_noint_ML, lme_full_ML);
disp(cmp);

% Save fitted model objects and outputs
% save(out_models, ...
%     'lme_full_ML', ...
%     'lme_noint_ML', ...
%     'cmp', ...
%     'a', ...
%     'f_full', ...
%     'f_noint', ...
%     'matFile', ...
%     '-v7.3');
% 
% fprintf('\nSaved model objects and statistical outputs to:\n%s\n', out_models);
% fprintf('\nLME analysis completed.\n');