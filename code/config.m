function paths = config()
% CONFIG  Local paths for the infant_energy pipeline.
%
%   Edit the fields below to point to your own copies of the data and
%   external toolboxes, then call:
%       paths = config();
%   from any script in this repository.

% ----------------------------------------------------------------------
% Data directories
% ----------------------------------------------------------------------
% Folder containing the cohort .mat files (dhcp_fetal.mat, dhcp_neo.mat,
% dhcp_fet_ce_alff.mat, dhcp_neo_ce_alff.mat, residuals.mat, etc.)
paths.data_dir = fullfile(pwd, '..', 'data');

% Folder containing BrainSpan transcriptomic data
% (expression_matrix.csv, columns_metadata.csv, rows_metadata.csv,
%  suppletable13.csv)
paths.brainspan_dir = fullfile(paths.data_dir, 'brainspan');

% Folder containing AHBA adult expression data (AHBA_adult_expression.csv)
paths.ahba_dir = fullfile(paths.data_dir, 'ahba');

% Folder containing the adult CMRglc atlas (adult_CMRglc.mat)
paths.cmrglc_dir = fullfile(paths.data_dir, 'cmrglc');

% Spin-test permutation file (infant_aal_spin_permutations_50k.mat)
paths.spin_perm_file = fullfile(paths.data_dir, ...
    'infant_aal_spin_permutations_50k.mat');

% ----------------------------------------------------------------------
% External toolboxes
% ----------------------------------------------------------------------
% Network control theory code (provides optim_fun.m).
% See: https://github.com/BassettLab/control_package
paths.nct_toolbox = '';   % e.g. '/path/to/control_package'

% Helper utilities (ALFF, spin_corr, scatter_corr_plot, ...).
% These live in code/helpers/ in this repository.
paths.helpers = fullfile(fileparts(mfilename('fullpath')), 'helpers');

% ----------------------------------------------------------------------
% Add everything to the MATLAB path
% ----------------------------------------------------------------------
if ~isempty(paths.nct_toolbox)
    addpath(genpath(paths.nct_toolbox));
end
addpath(genpath(paths.helpers));

end
