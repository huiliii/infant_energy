%% s04_transcriptomic_analysis.m
% Transcriptomic context for regions with postnatal increases in ALFF
% relative to control energy (Fig. 4).
%
% Two analyses:
%   (a) BrainSpan, 8-37 post-conception weeks: compare expression of
%       canonical neurodevelopmental gene sets (synapse development,
%       neuron migration, neuron differentiation, dendrite, axon,
%       myelination, cell proliferation) between regions of interest
%       and the rest of the brain.
%   (b) AHBA adult expression: per-gene t-tests over the left hemisphere
%       (45 left-hemisphere AAL parcels), exporting ranked gene lists
%       for downstream enrichment analysis.
%
% Inputs (configure paths in config.m):
%   <data_dir>/residuals.mat
%       ALFF_residuals_fetal, ALFF_residuals_preterm, ALFF_residuals_term
%       net_index (90 x 1; values 1-8 for AAL networks, 8 = subcortex)
%
%   <brainspan_dir>/
%       expression_matrix.csv
%       columns_metadata.csv
%       rows_metadata.csv
%       suppletable13.csv         (Kang et al. 2011 gene-set table)
%       trsp_brain_regions.mat    (aal_rns1, aal_rns2: AAL <-> BrainSpan
%                                   region mapping)
%
%   <ahba_dir>/AHBA_adult_expression.csv

clear; clc;

paths = config();

% ---------------------------------------------------------------------
%% Load ALFF residuals and define regions of interest
% ---------------------------------------------------------------------
load(fullfile(paths.data_dir, 'residuals.mat'), ...
    'ALFF_residuals_fetal', 'ALFF_residuals_preterm', ...
    'ALFF_residuals_term',  'net_index');

% Boolean maps of regions with positive ALFF residual per cohort.
brain_od_fetal   = (ALFF_residuals_fetal   > 0);
brain_od_preterm = (ALFF_residuals_preterm > 0);
brain_od_term    = (ALFF_residuals_term    > 0);  %#ok<NASGU>

% Regions of interest: subcortical network (net_index == 8). These are
% the regions with the largest postnatal increases in ALFF residuals.
brain_od_pos = (net_index == 8);

% Alternative definitions kept for reference:
%   brain_od_pos = (brain_od_fetal < brain_od_term);     % 0->1 transition
%   brain_od_pos = (ALFF_residuals_term > ALFF_residuals_fetal) ...
%                  & (ALFF_residuals_fetal < 0);

% ---------------------------------------------------------------------
%% BrainSpan analysis (developmental gene expression, 8-37 PCW)
% ---------------------------------------------------------------------
expression = readtable(fullfile(paths.brainspan_dir, 'expression_matrix.csv'));
sample     = readtable(fullfile(paths.brainspan_dir, 'columns_metadata.csv'));
rna        = readtable(fullfile(paths.brainspan_dir, 'rows_metadata.csv'));
gene_lists = readtable(fullfile(paths.brainspan_dir, 'suppletable13.csv'));

load(fullfile(paths.brainspan_dir, 'trsp_brain_regions.mat'), ...
    'aal_rns1', 'aal_rns2');

brain_regions = unique(sample.structure_name);

% Pick the gene set of interest. Columns in suppletable13:
%   Var2  - cell proliferation
%   Var6  - myelination
%   Var7  - axon development
%   Var8  - dendrite development
%   Var9  - neuron differentiation
%   Var10 - neuron migration
%   Var11 - synapse development
gene_set_name   = 'synapse development';
gene_set_column = gene_lists.Var11;

[~, gene_id] = intersect(rna.gene_symbol, gene_set_column);

% Drop the first column of expression_matrix.csv (row index, no header).
nd_related_gene = expression(gene_id, 2:525);

% Restrict to 8-37 post-conception weeks (first 237 samples).
num_sample      = 1:237;
expression_nd   = table2array(nd_related_gene(:, num_sample));
sample_regions  = sample.structure_name(num_sample);

% Average per BrainSpan region.
expression_nd_regional = zeros(size(expression_nd, 1), 26);
for r = 1:26
    idx = strcmp(sample_regions, brain_regions{r});
    expression_nd_regional(:, r) = mean(expression_nd(:, idx), 2);
end

% Map BrainSpan regions to the 90-node AAL parcellation
% (aal_rns1, aal_rns2 give two BrainSpan-region indices per AAL parcel).
brain_nd_rna_map = (expression_nd_regional(:, aal_rns1) + ...
                    expression_nd_regional(:, aal_rns2)) / 2;

% Compare gene expression in ROI vs rest of brain (paired across genes).
brain_change_nd_rna = brain_nd_rna_map(:, brain_od_pos == 1);
brain_keep_nd_rna   = brain_nd_rna_map(:, brain_od_pos == 0);

[h_bs, p_bs, ~, t_bs] = ttest(mean(brain_change_nd_rna, 2), ...
                              mean(brain_keep_nd_rna,   2));

fprintf('BrainSpan %s: t = %.3f, p = %.3g (h = %d)\n', ...
        gene_set_name, t_bs.tstat, p_bs, h_bs);

% ---------------------------------------------------------------------
%% Null-model permutation test (1000 random gene sets of matched size)
% ---------------------------------------------------------------------
run_null_model = false;   % flip to true to run (slow)

if run_null_model
    num_gene = size(nd_related_gene, 1);
    n_iter   = 1000;

    surrogate_gene_idx = randi([1 size(expression, 1)], num_gene, n_iter);
    p_permu = zeros(n_iter, 1);
    t_permu = zeros(n_iter, 1);

    for iter = 1:n_iter
        sur_gene  = expression(surrogate_gene_idx(:, iter), 2:525);
        sur_expr  = table2array(sur_gene(:, num_sample));
        sur_reg   = zeros(num_gene, 26);
        for r = 1:26
            idx = strcmp(sample_regions, brain_regions{r});
            sur_reg(:, r) = mean(sur_expr(:, idx), 2);
        end
        sur_map = (sur_reg(:, aal_rns1) + sur_reg(:, aal_rns2)) / 2;

        sur_change = sur_map(:, brain_od_pos == 1);
        sur_keep   = sur_map(:, brain_od_pos == 0);

        [~, p_permu(iter), ~, tt] = ttest(mean(sur_change, 2), ...
                                          mean(sur_keep,   2));
        t_permu(iter) = tt.tstat;
    end

    figure;
    histfit(p_permu, 100);
    hold on;
    xline(p_bs, 'r-', 'LineWidth', 3);
    xlabel('p'); ylabel('number of surrogate gene sets');
    title(gene_set_name); hold off;
end

% ---------------------------------------------------------------------
%% Validation with AHBA adult expression (left hemisphere only)
% ---------------------------------------------------------------------
AHBA_adult     = readtable(fullfile(paths.ahba_dir, 'AHBA_adult_expression.csv'));
AHBA_gene_list = AHBA_adult.Properties.VariableNames';
AHBA_gene_list = AHBA_gene_list(2:end);   % drop region-id column

% Restrict to left-hemisphere parcels (odd rows in AAL90).
AHBA_adult_left  = AHBA_adult(1:2:90, :);
AHBA_adult_left  = table2array(AHBA_adult_left(:, 2:end));
brain_od_pos_left = brain_od_pos(1:2:90);

brain_change_nd_rna = AHBA_adult_left(brain_od_pos_left == 1, :);
brain_keep_nd_rna   = AHBA_adult_left(brain_od_pos_left == 0, :);

num_gene_ahba = size(AHBA_adult_left, 2);
t_AHBA = zeros(num_gene_ahba, 1);
p_AHBA = zeros(num_gene_ahba, 1);

for g = 1:num_gene_ahba
    [~, p_AHBA(g), ~, tt] = ttest2(brain_change_nd_rna(:, g), ...
                                   brain_keep_nd_rna(:, g));
    t_AHBA(g) = tt.tstat;
end

% Significant genes (uncorrected); optional FDR / Bonferroni-Holm shown
% in commented lines.
%   p_AHBA_fdr = mafdr(p_AHBA, 'BHFDR', true);
%   p_AHBA_bh  = bonf_holm(p_AHBA);

sig_mask = p_AHBA < 0.05;
AHBA_sig = AHBA_gene_list(sig_mask);
writecell(AHBA_sig, fullfile(paths.data_dir, 'sig_gene_list.txt'));

% Split significant genes by direction and rank by |t|.
t_sig = t_AHBA;
t_sig(~sig_mask) = 0;

pos_mask = t_sig > 0;
neg_mask = t_sig < 0;

[~, ipos] = sort(t_sig(pos_mask), 'descend');
[~, ineg] = sort(t_sig(neg_mask), 'descend');

AHBA_pos_rank = AHBA_gene_list(pos_mask);  AHBA_pos_rank = AHBA_pos_rank(ipos);
AHBA_neg_rank = AHBA_gene_list(neg_mask);  AHBA_neg_rank = AHBA_neg_rank(ineg);

writecell(AHBA_pos_rank, fullfile(paths.data_dir, 'pos_gene_list.txt'));
writecell(AHBA_neg_rank, fullfile(paths.data_dir, 'neg_gene_list.txt'));

fprintf('AHBA significant genes (uncorrected p < 0.05): %d / %d\n', ...
        sum(sig_mask), num_gene_ahba);

% ---------------------------------------------------------------------
%% Optional: spin-test for the AHBA t-statistic
% ---------------------------------------------------------------------
% load(paths.spin_perm_file, 'perm_id');
% surrogate_brain_od_pos      = brain_od_pos(perm_id);
% surrogate_brain_od_pos_left = surrogate_brain_od_pos(1:2:90, :);
% sig_idx = find(sig_mask);
% n_spin  = size(perm_id, 2);
% t_AHBA_surr = zeros(numel(sig_idx), n_spin);
% p_AHBA_surr = zeros(numel(sig_idx), n_spin);
% for s = 1:n_spin
%     sc = AHBA_adult_left(surrogate_brain_od_pos_left(:, s) == 1, :);
%     sk = AHBA_adult_left(surrogate_brain_od_pos_left(:, s) == 0, :);
%     for g = 1:numel(sig_idx)
%         [~, p_AHBA_surr(g, s), ~, tt] = ttest2(sc(:, sig_idx(g)), ...
%                                                sk(:, sig_idx(g)));
%         t_AHBA_surr(g, s) = tt.tstat;
%     end
% end
