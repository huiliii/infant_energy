function [CE, ALFF] = compute_subject_ce_alff(A, tc, TR, num_timepoint, num_node)
% COMPUTE_SUBJECT_CE_ALFF  Per-subject control energy and ALFF.
%
%   [CE, ALFF] = compute_subject_ce_alff(A, tc, TR, num_timepoint, num_node)
%
%   Inputs
%     A             : num_node x num_node structural connectivity matrix
%     tc            : num_timepoint x num_node ROI time series
%     TR            : repetition time of the fMRI sequence (seconds)
%     num_timepoint : number of fMRI volumes
%     num_node      : number of parcels (90 for infant AAL)
%
%   Outputs
%     CE   : num_node x 1 control energy per region (optimal control,
%            baseline -> all-ones target state)
%     ALFF : num_node x 1 amplitude of low-frequency fluctuations
%
%   Requires `optim_fun` (network control theory toolbox) and `cal_alff`
%   (helper). See README for sources.

% ----- Control energy -------------------------------------------------
% Normalise A by its largest eigenvalue (+1) and subtract identity so
% the continuous-time system x'(t) = A x(t) + B u(t) is stable.
A_norm = A ./ (eigs(A, 1) + 1) - eye(size(A, 1));

base_state   = zeros(num_node, 1);
target_state = ones(num_node, 1);

[~, U_opt, ~] = optim_fun(A_norm, 1, base_state, target_state, 1);
CE = trapz(U_opt.^2)';

% ----- ALFF -----------------------------------------------------------
[ALFF, ~] = cal_alff(tc, TR, num_timepoint, num_node);

end
