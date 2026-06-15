%% s02_fig2_ce_alff_coupling.m
% Reproduces the analyses summarised in Figure 2:
%   - whole-brain CE ~ age within each cohort
%   - spatial similarity of group-averaged CE and ALFF maps across cohorts
%   - group-averaged CE-ALFF spatial correlation per cohort
%   - group differences in individual-level CE-ALFF correlations
%   - r(CE, ALFF) ~ age, with and without covariates
%
% Requires: s01_compute_ce_alff.m to have been run, plus
%           - perm_id from infant_aal_spin_permutations_50k.mat
%           - tsnr_z_fet, tsnr_z_neo (fMRI motion proxy)
%           - pct_slice (fetal dMRI motion), shard_outlier (neo dMRI motion)
%
% Helper functions used: spin_corr, scatter_corr_plot.

clear; clc;

paths    = config();
data_dir = paths.data_dir;

% ---------------------------------------------------------------------
%% Load CE, ALFF, ages, motion, spin permutations
% ---------------------------------------------------------------------
load(fullfile(data_dir, 'dhcp_fet_ce_alff.mat'), 'CE_opt_wb_fet', 'ALFF_fet');
load(fullfile(data_dir, 'dhcp_neo_ce_alff.mat'), 'CE_opt_wb_neo', 'ALFF_neo');

load(fullfile(data_dir, 'dhcp_fetal.mat'), 'age', 'tsnr_z_fet', 'pct_slice');
load(fullfile(data_dir, 'dhcp_neo.mat'),   'pma', 'index_term', ...
                                            'tsnr_z_neo', 'shard_outlier');

load(fullfile(data_dir, 'group_corr.mat'), 'r_fet', 'r_neo');

load(paths.spin_perm_file, 'perm_id');

% ---------------------------------------------------------------------
%% Whole-brain control energy vs age, within each cohort
% ---------------------------------------------------------------------
[r_age_fet,     p_age_fet]     = corr(mean(CE_opt_wb_fet)', age);
[r_age_preterm, p_age_preterm] = corr(mean(CE_opt_wb_neo(:, ~index_term))', ...
                                       pma(~index_term));
[r_age_term,    p_age_term]    = corr(mean(CE_opt_wb_neo(:,  index_term))', ...
                                       pma( index_term));

figure;
subplot(1, 3, 1);
scatter_corr_plot(age, mean(CE_opt_wb_fet)');
xlabel('fetal age / week'); ylabel('control energy');

subplot(1, 3, 2);
scatter_corr_plot(pma(~index_term), mean(CE_opt_wb_neo(:, ~index_term))');
xlabel('preterm age / week'); ylabel('control energy');

subplot(1, 3, 3);
scatter_corr_plot(pma(index_term), mean(CE_opt_wb_neo(:, index_term))');
xlabel('term age / week'); ylabel('control energy');

% ---------------------------------------------------------------------
%% Spatial similarity of group-averaged CE maps across cohorts
% ---------------------------------------------------------------------
mean_CE_fet     = mean(CE_opt_wb_fet, 2);
mean_CE_preterm = mean(CE_opt_wb_neo(:, ~index_term), 2);
mean_CE_term    = mean(CE_opt_wb_neo(:,  index_term), 2);

[r_ce_fp, p_ce_fp, p_ce_fp_spin] = spin_corr(mean_CE_fet,     mean_CE_preterm, perm_id);
[r_ce_ft, p_ce_ft, p_ce_ft_spin] = spin_corr(mean_CE_fet,     mean_CE_term,    perm_id);
[r_ce_pt, p_ce_pt, p_ce_pt_spin] = spin_corr(mean_CE_preterm, mean_CE_term,    perm_id);

% ---------------------------------------------------------------------
%% Spatial similarity of group-averaged ALFF maps across cohorts
% ---------------------------------------------------------------------
mean_ALFF_fet     = mean(ALFF_fet, 2);
mean_ALFF_preterm = mean(ALFF_neo(:, ~index_term), 2);
mean_ALFF_term    = mean(ALFF_neo(:,  index_term), 2);

[r_al_fp, p_al_fp, p_al_fp_spin] = spin_corr(mean_ALFF_fet,     mean_ALFF_preterm, perm_id);
[r_al_ft, p_al_ft, p_al_ft_spin] = spin_corr(mean_ALFF_fet,     mean_ALFF_term,    perm_id);
[r_al_pt, p_al_pt, p_al_pt_spin] = spin_corr(mean_ALFF_preterm, mean_ALFF_term,    perm_id);

% ---------------------------------------------------------------------
%% Group-averaged CE vs ALFF spatial correlations (Fig. 2a-c)
% ---------------------------------------------------------------------
[r_ca_fet,     p_ca_fet,     p_ca_fet_spin]     = spin_corr(mean_CE_fet,     mean_ALFF_fet,     perm_id);
[r_ca_preterm, p_ca_preterm, p_ca_preterm_spin] = spin_corr(mean_CE_preterm, mean_ALFF_preterm, perm_id);
[r_ca_term,    p_ca_term,    p_ca_term_spin]    = spin_corr(mean_CE_term,    mean_ALFF_term,    perm_id);

% ---------------------------------------------------------------------
%% Group differences in individual-level CE-ALFF correlations (Fig. 2d)
% ---------------------------------------------------------------------
[~, p_ft, ~, t_ft] = ttest2(r_fet,             r_neo( index_term));
[~, p_pt, ~, t_pt] = ttest2(r_neo(~index_term), r_neo( index_term));
[~, p_fp, ~, t_fp] = ttest2(r_fet,             r_neo(~index_term));

% ---------------------------------------------------------------------
%% r(CE, ALFF) ~ age across the whole perinatal period (Fig. 2e)
% ---------------------------------------------------------------------
r_all   = [r_fet;       r_neo];
age_all = [age;         pma];

[r_age, p_age] = corr(r_all, age_all);

% Partial correlation controlling for fMRI tSNR and dMRI motion.
tsnr_z     = [tsnr_z_fet; tsnr_z_neo];
dti_motion = [pct_slice;  shard_outlier];
covariates = [tsnr_z, dti_motion];

[r_age_partial, p_age_partial] = partialcorr(r_all, age_all, covariates);

fprintf('\nCE-ALFF coupling vs age:\n');
fprintf('  raw      : r = %.3f, p = %.3g\n', r_age,         p_age);
fprintf('  partial  : r = %.3f, p = %.3g  (controls: tSNR, dMRI motion)\n', ...
        r_age_partial, p_age_partial);
