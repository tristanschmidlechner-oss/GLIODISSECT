% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function test_normative_model
rng(20260909);
G = linspace(-0.09, -0.04, 20)';
beta0 = -0.02; beta1 = 0.5;
noise = 0.001 * randn(2,2,1,20);
maps = beta0 + beta1 * reshape(G,1,1,1,[]) + noise;
anatomy = true(2,2,1);
subjectMasks = true(size(maps));
model = fit_g_adjusted_model(maps, G, anatomy, subjectMasks, 0.8);
assert(all(abs(model.beta0(:) - beta0) < 0.01));
assert(all(abs(model.beta1(:) - beta1) < 0.2));
assert(all(model.r_squared(:) > 0));

observed = beta0 + beta1 * G(1) + zeros(2,2,1);
result = apply_g_adjusted_model(observed, G(1), anatomy, model, 2);
assert(all(isfinite(result.expected(:))));
assert(max(abs(result.difference(:))) < 0.01);
end
