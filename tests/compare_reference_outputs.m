% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function metrics = compare_reference_outputs(generated, reference)
%COMPARE_REFERENCE_OUTPUTS Compare generated and archived map structures.
%
% Load two MAT files containing structures with identical field names. The
% function reports finite-voxel Pearson correlation, RMSE and median
% absolute difference for each numeric map.

generatedData = load(generated); referenceData = load(reference);
generatedMaps = select_maps(generatedData); referenceMaps = select_maps(referenceData);
names = intersect(fieldnames(generatedMaps), fieldnames(referenceMaps), 'stable');
rows = cell(0,5);
for index = 1:numel(names)
    name = names{index}; A = generatedMaps.(name); B = referenceMaps.(name);
    if ~isnumeric(A) || ~isequal(size(A), size(B)); continue; end
    keep = isfinite(A) & isfinite(B);
    if nnz(keep) < 2; continue; end
    a = double(A(keep)); b = double(B(keep));
    R = corrcoef(a,b);
    rows(end+1,:) = {name, nnz(keep), R(1,2), sqrt(mean((a-b).^2)), median(abs(a-b))}; %#ok<AGROW>
end
metrics = cell2table(rows, 'VariableNames', ...
    {'map','common_valid_voxels','pearson_r','rmse','median_absolute_difference'});
end

function maps = select_maps(data)
if isfield(data, 'maps'); maps = data.maps; else; maps = data; end
end
