% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function mask = build_healthy_mask(greyMatter, whiteMatter, CSF, input, parameters)
%BUILD_HEALTHY_MASK Construct the healthy-control EPI analysis mask.

validate_same_size(greyMatter, whiteMatter, CSF);
threshold = parameters.tissue_probability_threshold;
if isfield(input, 'tissueProbabilityThreshold')
    threshold = input.tissueProbabilityThreshold;
end
mask = isfinite(greyMatter) & isfinite(whiteMatter) & isfinite(CSF) & ...
    (greyMatter + whiteMatter + CSF >= threshold);
end

function validate_same_size(varargin)
referenceSize = size(varargin{1});
if any(~cellfun(@(x) isequal(size(x), referenceSize), varargin))
    error('HypoxiaBOLD:MaskSizeMismatch', 'All tissue maps must have identical dimensions.');
end
end
