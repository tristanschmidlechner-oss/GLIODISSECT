% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function test_filter_shape
rng(20260909);
raw = 100 + randn(2,2,1,160);
mask = true(2,2,1);
filtered = filter_bold_timeseries(raw, mask, 1.8, 12);
assert(isequal(size(filtered), size(raw)));
assert(all(isfinite(filtered(:))));
end
