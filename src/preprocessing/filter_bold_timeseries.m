% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function filteredBOLD = filter_bold_timeseries(rawBOLD, analysisMask, TR, windowVolumes)
%FILTER_BOLD_TIMESERIES Apply the manuscript temporal preprocessing.
%
% The frequency-domain filter is followed by robust locally weighted
% smoothing over a fixed number of volumes. The default publication
% configuration is 12 volumes (21.6 s at TR = 1.8 s). This implementation
% uses base-MATLAB smoothdata and therefore does not require Curve Fitting
% Toolbox.

if nargin < 3 || isempty(TR); TR = 1.8; end
if nargin < 4 || isempty(windowVolumes); windowVolumes = 12; end
validateattributes(rawBOLD, {'numeric'}, {'real', 'nonempty'});
spatialSize = [size(rawBOLD, 1), size(rawBOLD, 2), size(rawBOLD, 3)];
validateattributes(analysisMask, {'logical','numeric'}, {'size', spatialSize});
validateattributes(TR, {'numeric'}, {'scalar','positive','finite'});
validateattributes(windowVolumes, {'numeric'}, {'scalar','integer','>=',3});

numberOfScans = size(rawBOLD, 4);
voxelSeries = reshape(rawBOLD, [], numberOfScans);
maskVector = logical(analysisMask(:));
filteredSeries = nan(size(voxelSeries));

for voxelIndex = find(maskVector)'
    series = double(voxelSeries(voxelIndex, :));
    if any(~isfinite(series)); continue; end
    detrended = detrend(series(:)) + mean(series);
    frequencyFiltered = low_highPassFilter_g(detrended, ...
        [0.125, 0, TR, numberOfScans], 'smooth box');
    smoothed = smoothdata(frequencyFiltered(:), 'rloess', windowVolumes);
    filteredSeries(voxelIndex, :) = smoothed(:)';
end

filteredBOLD = reshape(filteredSeries, size(rawBOLD));
end
