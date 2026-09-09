% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function result = apply_g_adjusted_model(observedMap, patientG, patientMask, model, zThreshold)
%APPLY_G_ADJUSTED_MODEL Create expected, difference and deviation maps.

if nargin < 5 || isempty(zThreshold); zThreshold = 2; end
required = {'beta0','beta1','residual_sd_raw','residual_sd_floor','eligible_mask'};
missing = required(~isfield(model, required));
if ~isempty(missing)
    error('HypoxiaBOLD:IncompleteModel', 'Missing model fields: %s', strjoin(missing, ', '));
end
validateattributes(patientMask, {'logical','numeric'}, {'size', size(observedMap)});

valid = logical(patientMask) & logical(model.eligible_mask) & ...
    isfinite(observedMap) & isfinite(model.beta0) & isfinite(model.beta1) & ...
    isfinite(model.residual_sd_raw);
expected = nan(size(observedMap)); difference = expected; z = expected;
expected(valid) = model.beta0(valid) + model.beta1(valid) .* patientG;
difference(valid) = observedMap(valid) - expected(valid);
denominator = max(model.residual_sd_raw(valid), model.residual_sd_floor);
z(valid) = difference(valid) ./ denominator;

signedAbnormal = zeros(size(observedMap), 'int8');
signedAbnormal(valid & z < -zThreshold) = -1;
signedAbnormal(valid & z > zThreshold) = 1;

result = struct('expected', expected, 'difference', difference, 'z', z, ...
    'signed_abnormal', signedAbnormal, 'valid_mask', valid);
end
