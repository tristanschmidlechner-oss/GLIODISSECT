% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function test_masks
gm = zeros(2,2,1); wm = gm; csf = gm;
gm(1) = 0.5; wm(1) = 0.3;
gm(2) = 0.2; wm(2) = 0.2;
lesion = false(size(gm)); lesion(2) = true;
parameters = default_parameters();
input = struct();
healthy = build_healthy_mask(gm, wm, csf, input, parameters);
patient = build_patient_mask(gm, wm, csf, lesion, input, parameters);
assert(healthy(1) && ~healthy(2));
assert(patient(1) && patient(2));
end
