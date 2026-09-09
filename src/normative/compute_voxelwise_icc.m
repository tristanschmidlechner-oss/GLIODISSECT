% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function result = compute_voxelwise_icc(measurement1, measurement2, subjectMasks, minimumCoverage)
%COMPUTE_VOXELWISE_ICC Compute ICC(2,1) and ICC(3,1) for paired maps.
%
% measurement1 and measurement2 are X-by-Y-by-Z-by-subject arrays.
% ICC(2,1) measures absolute agreement; ICC(3,1) measures consistency.

if nargin < 4 || isempty(minimumCoverage); minimumCoverage = 0.80; end
if ~isequal(size(measurement1), size(measurement2), size(subjectMasks))
    error('HypoxiaBOLD:ICCSizeMismatch', 'Measurements and subject masks must have identical dimensions.');
end
numberOfSubjects = size(measurement1,4);
minimumN = ceil(minimumCoverage * numberOfSubjects);
spatialSize = [size(measurement1,1), size(measurement1,2), size(measurement1,3)];
A = reshape(double(measurement1), [], numberOfSubjects);
B = reshape(double(measurement2), [], numberOfSubjects);
M = reshape(logical(subjectMasks), [], numberOfSubjects);
icc21 = nan(size(A,1),1); icc31 = icc21; pearsonR = icc21; sampleSize = zeros(size(A,1),1);

for voxelIndex = 1:size(A,1)
    keep = M(voxelIndex,:) & isfinite(A(voxelIndex,:)) & isfinite(B(voxelIndex,:));
    n = sum(keep);
    if n < max(minimumN,3); continue; end
    Y = [A(voxelIndex,keep)', B(voxelIndex,keep)'];
    grandMean = mean(Y,'all');
    rowMeans = mean(Y,2); columnMeans = mean(Y,1);
    k = 2;
    SSR = k * sum((rowMeans - grandMean).^2);
    SSC = n * sum((columnMeans - grandMean).^2);
    residual = Y - rowMeans - columnMeans + grandMean;
    SSE = sum(residual.^2,'all');
    MSR = SSR / (n - 1);
    MSC = SSC / (k - 1);
    MSE = SSE / ((n - 1) * (k - 1));
    denominator21 = MSR + (k - 1) * MSE + k * (MSC - MSE) / n;
    denominator31 = MSR + (k - 1) * MSE;
    if denominator21 ~= 0; icc21(voxelIndex) = (MSR - MSE) / denominator21; end
    if denominator31 ~= 0; icc31(voxelIndex) = (MSR - MSE) / denominator31; end
    R = corrcoef(Y(:,1),Y(:,2)); pearsonR(voxelIndex) = R(1,2);
    sampleSize(voxelIndex) = n;
end

result = struct('ICC_2_1',reshape(icc21,spatialSize), ...
    'ICC_3_1',reshape(icc31,spatialSize), ...
    'pearson_r',reshape(pearsonR,spatialSize), ...
    'sample_size',reshape(sampleSize,spatialSize), ...
    'minimum_coverage_fraction',minimumCoverage);
end
