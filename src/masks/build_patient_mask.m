% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function mask = build_patient_mask(greyMatter, whiteMatter, CSF, lesionMask, input, parameters)
%BUILD_PATIENT_MASK Construct the patient EPI mask from segmentation output.
%
% Patient lesions are retained even where tissue segmentation probability
% is reduced. This adapter is intentionally distinct from the healthy mask.

referenceSize = size(greyMatter);
inputs = {whiteMatter, CSF, lesionMask};
if any(~cellfun(@(x) isequal(size(x), referenceSize), inputs))
    error('HypoxiaBOLD:MaskSizeMismatch', 'Tissue and lesion maps must have identical dimensions.');
end
threshold = parameters.tissue_probability_threshold;
if isfield(input, 'tissueProbabilityThreshold')
    threshold = input.tissueProbabilityThreshold;
end
tissueMask = isfinite(greyMatter) & isfinite(whiteMatter) & isfinite(CSF) & ...
    (greyMatter + whiteMatter + CSF >= threshold);
mask = tissueMask | (isfinite(lesionMask) & lesionMask > 0);
end
