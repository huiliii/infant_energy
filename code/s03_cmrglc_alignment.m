%% s03_fig5_cmrglc_alignment.m
% Reproduces Figure 5: spatial similarity between individual control
% energy maps and the adult cerebral metabolic rate of glucose (CMRglc)
% map increases across the perinatal period.
%
% The adult CMRglc map comes from Shokri-Kojori et al. 2019 (28 adults),
% registered to the infant AAL90 parcellation.
%
% Inputs (in <data_dir>/):
%   dhcp_fet_ce_alff.mat   - CE_opt_wb_fet
%   dhcp_neo_ce_alff.mat   - CE_opt_wb_neo
%   dhcp_fetal.mat         - age
%   dhcp_neo.mat           - pma
%   <cmrglc_dir>/adult_CMRglc.mat - adult_CMRglc_aal (90 x 1)

clear; clc;

paths = config();

load(fullfile(paths.data_dir,   'dhcp_fet_ce_alff.mat'), 'CE_opt_wb_fet');
load(fullfile(paths.data_dir,   'dhcp_neo_ce_alff.mat'), 'CE_opt_wb_neo');
load(fullfile(paths.data_dir,   'dhcp_fetal.mat'),       'age');
load(fullfile(paths.data_dir,   'dhcp_neo.mat'),         'pma');
load(fullfile(paths.cmrglc_dir, 'adult_CMRglc.mat'),     'adult_CMRglc_aal');

% Concatenate fetal + neonatal cohorts in subject order.
CE_opt_wb = cat(2, CE_opt_wb_fet, CE_opt_wb_neo);
age_all   = [age; pma];

num_sub = numel(age_all);

% Per-subject spatial correlation of CE map with the adult CMRglc map.
rCC = zeros(num_sub, 1);
pCC = zeros(num_sub, 1);
for s = 1:num_sub
    [rCC(s), pCC(s)] = corr(CE_opt_wb(:, s), adult_CMRglc_aal);
end

% Does this per-subject CE-CMRglc spatial correlation grow with age?
[r_age, p_age] = corr(rCC, age_all);

fprintf('r(CE,CMRglc) vs age: r = %.3f, p = %.3g\n', r_age, p_age);
