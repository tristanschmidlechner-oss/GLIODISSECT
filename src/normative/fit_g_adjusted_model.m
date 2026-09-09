% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function model = fit_g_adjusted_model(controlMaps, G, anatomicalMask, subjectMasks, minimumCoverage)
%FIT_G_ADJUSTED_MODEL Fit the corrected voxelwise normative OLS model.
%
% controlMaps has dimensions X-by-Y-by-Z-by-subject. Anatomical eligibility
% and subject-level masks are kept separate from data coverage. R-squared is
% ordinary in-sample OLS R-squared and is not used for patient prediction.

if nargin < 5 || isempty(minimumCoverage); minimumCoverage = 0.80; end
numberOfSubjects = size(controlMaps, 4);
G = double(G(:));
validateattributes(G, {'numeric'}, {'vector','numel',numberOfSubjects,'finite'});
spatialSize = [size(controlMaps, 1), size(controlMaps, 2), size(controlMaps, 3)];
validateattributes(anatomicalMask, {'logical','numeric'}, {'size', spatialSize});
validateattributes(subjectMasks, {'logical','numeric'}, {'size', size(controlMaps)});

Y = reshape(double(controlMaps), [], numberOfSubjects);
M = reshape(logical(subjectMasks), [], numberOfSubjects);
anatomy = logical(anatomicalMask(:));
valid = M & isfinite(Y);
minimumN = ceil(minimumCoverage * numberOfSubjects);
eligible = anatomy & sum(valid, 2) >= minimumN;

beta0 = nan(size(Y,1),1); beta1 = beta0; residualSD = beta0; rSquared = beta0;
for voxelIndex = find(eligible)'
    keep = valid(voxelIndex,:)' & isfinite(G);
    y = Y(voxelIndex,keep)';
    X = [ones(sum(keep),1), G(keep)];
    beta = X \ y;
    fitted = X * beta;
    residual = y - fitted;
    SSE = sum(residual.^2);
    SST = sum((y - mean(y)).^2);
    beta0(voxelIndex) = beta(1);
    beta1(voxelIndex) = beta(2);
    residualSD(voxelIndex) = sqrt(SSE / (numel(y) - 2));
    if SST > 0
        rSquared(voxelIndex) = 1 - SSE / SST;
    end
end

positiveSD = residualSD(eligible & isfinite(residualSD) & residualSD > 0);
if isempty(positiveSD)
    error('HypoxiaBOLD:NoResidualVariance', 'No positive residual standard deviations were estimated.');
end
floorValue = prctile(positiveSD, 5);

model = struct();
model.beta0 = reshape(beta0, spatialSize);
model.beta1 = reshape(beta1, spatialSize);
model.residual_sd_raw = reshape(residualSD, spatialSize);
model.residual_sd_floor = floorValue;
model.r_squared = reshape(rSquared, spatialSize);
model.eligible_mask = reshape(eligible, spatialSize);
model.minimum_coverage_fraction = minimumCoverage;
model.minimum_subject_count = minimumN;
model.number_of_controls = numberOfSubjects;
end
