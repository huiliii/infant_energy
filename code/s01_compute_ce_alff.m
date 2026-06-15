%% 01_compute_ce_alff.m
% Compute regional control energy and ALFF for each subject in the
% fetal and neonatal dHCP cohorts, and the within-subject spatial
% correlation between them (Fig. 1, Fig. 2).
%
% Inputs (in <data_dir>/):
%   dhcp_fetal.mat   - struct fields: mat (90x90xN), tcs (350x90xN), age
%   dhcp_neo.mat     - struct fields: mat (90x90xN), tcs (2300x90xN), pma,
%                      index_term (logical, 1 for term)
%
% Outputs (saved to <data_dir>/):
%   dhcp_fet_ce_alff.mat  - CE_opt_wb_fet, ALFF_fet
%   dhcp_neo_ce_alff.mat  - CE_opt_wb_neo, ALFF_neo
%   group_corr.mat        - r_fet, p_fet, r_neo, p_neo (per-subject corr),
%                           FDR-corrected versions, and group-level corr

clear; clc;

paths = config();
data_dir = paths.data_dir;

% Whether to recompute CE+ALFF from the connectivity matrices and time
% series, or load pre-computed values. Set to false to skip the heavy
% optimisation loop.
do_compute = true;

% ---------------------------------------------------------------------
%% Fetal cohort
% ---------------------------------------------------------------------
load(fullfile(data_dir, 'dhcp_fetal.mat'), 'mat', 'tcs');

TR_fet            = 2.2;     % seconds
num_timepoint_fet = 350;

if do_compute
    [CE_opt_wb_fet, ALFF_fet] = compute_cohort_ce_alff( ...
        mat, tcs, TR_fet, num_timepoint_fet);
    save(fullfile(data_dir, 'dhcp_fet_ce_alff.mat'), ...
        'CE_opt_wb_fet', 'ALFF_fet');
else
    load(fullfile(data_dir, 'dhcp_fet_ce_alff.mat'), ...
        'CE_opt_wb_fet', 'ALFF_fet');
end

num_fet = size(CE_opt_wb_fet, 2);

% ---------------------------------------------------------------------
%% Neonatal cohort (preterm + term)
% ---------------------------------------------------------------------
load(fullfile(data_dir, 'dhcp_neo.mat'), 'mat', 'tcs');

TR_neo            = 0.392;
num_timepoint_neo = 2300;

if do_compute
    [CE_opt_wb_neo, ALFF_neo] = compute_cohort_ce_alff( ...
        mat, tcs, TR_neo, num_timepoint_neo);
    save(fullfile(data_dir, 'dhcp_neo_ce_alff.mat'), ...
        'CE_opt_wb_neo', 'ALFF_neo');
else
    load(fullfile(data_dir, 'dhcp_neo_ce_alff.mat'), ...
        'CE_opt_wb_neo', 'ALFF_neo');
end

num_neo = size(CE_opt_wb_neo, 2);

% ---------------------------------------------------------------------
%% Group-level and individual-level CE-ALFF spatial correlations
% ---------------------------------------------------------------------
% Group-level: correlate the group-averaged CE and ALFF maps.
[r_fet_group, p_fet_group] = corr(mean(ALFF_fet,       2), ...
                                  mean(CE_opt_wb_fet,  2));
[r_neo_group, p_neo_group] = corr(mean(ALFF_neo,       2), ...
                                  mean(CE_opt_wb_neo,  2));

% Individual-level: per-subject spatial correlation across the 90 ROIs.
r_fet = zeros(num_fet, 1);  p_fet = zeros(num_fet, 1);
for s = 1:num_fet
    [r_fet(s), p_fet(s)] = corr(ALFF_fet(:, s), CE_opt_wb_fet(:, s));
end

r_neo = zeros(num_neo, 1);  p_neo = zeros(num_neo, 1);
for s = 1:num_neo
    [r_neo(s), p_neo(s)] = corr(ALFF_neo(:, s), CE_opt_wb_neo(:, s));
end

% Benjamini-Hochberg FDR correction within each cohort.
p_fetFDR = mafdr(p_fet, 'BHFDR', true);
p_neoFDR = mafdr(p_neo, 'BHFDR', true);

r_fetFDR                = r_fet;
r_fetFDR(p_fetFDR > 0.05) = 0;
r_neoFDR                = r_neo;
r_neoFDR(p_neoFDR > 0.05) = 0;

save(fullfile(data_dir, 'group_corr.mat'), ...
    'r_fet_group', 'p_fet_group', 'r_neo_group', 'p_neo_group', ...
    'r_fet', 'p_fet', 'r_fetFDR', 'p_fetFDR', ...
    'r_neo', 'p_neo', 'r_neoFDR', 'p_neoFDR');

fprintf('Fetal group-level   r(CE,ALFF) = %.3f, p = %.3g\n', ...
        r_fet_group, p_fet_group);
fprintf('Neonatal group-level r(CE,ALFF) = %.3f, p = %.3g\n', ...
        r_neo_group, p_neo_group);
fprintf('Fetuses with FDR-significant individual coupling: %d / %d\n', ...
        sum(p_fetFDR < 0.05), num_fet);
fprintf('Neonates with FDR-significant individual coupling: %d / %d\n', ...
        sum(p_neoFDR < 0.05), num_neo);
