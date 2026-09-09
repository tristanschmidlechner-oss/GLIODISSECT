% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function test_icc
rng(20260909);
n = 20;
subjectEffect = reshape(linspace(-1,1,n),1,1,1,[]);
measurement1 = subjectEffect + 0.01 * randn(2,2,1,n);
measurement2 = subjectEffect + 0.01 * randn(2,2,1,n);
result = compute_voxelwise_icc(measurement1,measurement2,true(size(measurement1)),0.8);
assert(all(result.ICC_2_1(:) > 0.95));
assert(all(result.ICC_3_1(:) > 0.95));
end
