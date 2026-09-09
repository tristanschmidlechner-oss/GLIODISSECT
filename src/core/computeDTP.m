% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [DTPmap, stepAmplitude, index10map, index90map, delayCorrection, Meanbaseline] = ...
    computeDTP(BOLD4D, brainlag, DurationRA, lagDiff2, Rowscans)
% computeDTP - Computes 10-90% rise time (DTP) and related maps for each voxel
%
% INPUTS:
%   BOLD4D       - 4D filtered fMRI data (e.g., bBOLD), size: X x Y x Z x T
%   brainlag     - 3D lag map (e.g., brainlagDiffO2), same X x Y x Z
%   DurationRA   - Vector [baselineDuration, stepDuration, maxDuration]
%   lagDiff2     - Offset value used in baseline shifting
%   Rowscans     - Time vector of scan TRs (e.g., [1, 2, 3, ..., N])
%
% OUTPUTS:
%   DTPmap           - Map of rise time from 10% to 90% (index90 - index10)
%   stepAmplitude    - BOLD step change amplitude
%   index10map       - Map of 10% response time indices
%   index90map       - Map of 90% response time indices
%   delayCorrection  - Optional delay correction map

% Get matrix size
[MtxX, MtxY, MtxZ, ~] = size(BOLD4D);

% Preallocate outputs
DTPmap = nan(MtxX, MtxY, MtxZ);
stepAmplitude = nan(MtxX, MtxY, MtxZ);
index10map = nan(MtxX, MtxY, MtxZ);
index90map = nan(MtxX, MtxY, MtxZ);
delayCorrection = nan(MtxX, MtxY, MtxZ);
Meanbaseline = nan(MtxX, MtxY, MtxZ);

% Extract durations
Durationbaseline = DurationRA(1);
Durationstepchange = DurationRA(2);
maxDuration = DurationRA(3);

% Loop over all voxels
for ix = 1:MtxX
    for iy = 1:MtxY
        for iz = 1:MtxZ
            shift = brainlag(ix, iy, iz);

            % Skip invalid voxels
            if isnan(shift) || shift > maxDuration
                continue
            end

            % Adjust baseline for voxel-specific lag
            Durationbaseline1 = Durationbaseline - lagDiff2 + shift;
            fData = squeeze(BOLD4D(ix, iy, iz, :));

            try
                [finalindex90, finalindex10, finalmeanbaseline, finalstep] = ...
                    complrange10_90c(Durationbaseline1, Durationstepchange, Rowscans, fData, shift);

                % Compute delay correction
                baselineShifted = fix(Durationbaseline1) - finalindex10;
                if baselineShifted >= shift
                    delayCorr = 0;
                else
                    delayCorr = shift - baselineShifted;
                end

                % Store outputs
                DTPmap(ix, iy, iz) = finalindex90 - finalindex10;
                stepAmplitude(ix, iy, iz) = finalstep;
                index10map(ix, iy, iz) = finalindex10;
                index90map(ix, iy, iz) = finalindex90;
                delayCorrection(ix, iy, iz) = delayCorr;
                Meanbaseline(ix, iy, iz) = finalmeanbaseline;

            catch
                % Leave NaN on error
                continue
            end
        end
    end
end
end
