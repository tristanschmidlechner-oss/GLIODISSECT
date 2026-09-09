% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function maps = run_healthy_control(input)
%RUN_HEALTHY_CONTROL Run the shared analysis with the healthy-control mask.
%
% The input structure must contain bBOLD, brainLag, durationRA,
% initialGlobalLag, rowScans, shiftedO2, greyMatter, whiteMatter and CSF.

required = {'bBOLD','brainLag','durationRA','initialGlobalLag','rowScans', ...
    'shiftedO2','greyMatter','whiteMatter','CSF'};
assert_fields(input, required);
parameters = default_parameters();
mask = build_healthy_mask(input.greyMatter, input.whiteMatter, input.CSF, input, parameters);
maps = run_hypoxia_bold_core(input.bBOLD, input.brainLag, input.durationRA, ...
    input.initialGlobalLag, input.rowScans, input.shiftedO2, mask);
maps.cohort = 'healthy_control';
end

function assert_fields(input, required)
missing = required(~isfield(input, required));
if ~isempty(missing)
    error('HypoxiaBOLD:MissingInput', 'Missing input fields: %s', strjoin(missing, ', '));
end
end
