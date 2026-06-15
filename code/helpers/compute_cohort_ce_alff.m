function [CE_all, ALFF_all] = compute_cohort_ce_alff(mat, tcs, TR, num_timepoint)
% COMPUTE_COHORT_CE_ALFF  Loop compute_subject_ce_alff over a cohort.
%
%   [CE_all, ALFF_all] = compute_cohort_ce_alff(mat, tcs, TR, num_timepoint)
%
%   Inputs
%     mat           : num_node x num_node x num_subject connectivity stack
%     tcs           : num_timepoint x num_node x num_subject ROI time series
%     TR            : repetition time (seconds)
%     num_timepoint : number of fMRI volumes
%
%   Outputs
%     CE_all   : num_node x num_subject control energy
%     ALFF_all : num_node x num_subject ALFF

[num_node, ~, num_subject] = size(mat);

CE_all   = zeros(num_node, num_subject);
ALFF_all = zeros(num_node, num_subject);

for s = 1:num_subject
    [CE_all(:, s), ALFF_all(:, s)] = compute_subject_ce_alff( ...
        mat(:, :, s), tcs(:, :, s), TR, num_timepoint, num_node);
end

end
