% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [DTRmap, index90_DTR_map, index10_DTR_map, Meanbaseline2map] = ...
    computeDTR(BOLD4D, brainlag, DurationRA, lagDiff2, Rowscans, index90_DTP, stepAmplitude)
% computeDTR  Voxelwise Delay-To-Return (DTR) for the O2 BOLD protocol
%
% DTR quantifies how quickly the BOLD signal returns from plateau back
% toward the post-step baseline after the stimulus step ends.
%
% Key indices:
%   index90_DTR  = last TR where signal is still >= 90% of plateau amplitude
%                  --> PLATEAU END (return has not yet started)
%   index10_DTR  = first TR where signal has returned to within 10% of
%                  post-step baseline level
%                  --> RETURN NEARLY COMPLETE
%
% Together with computeDTP outputs, the per-voxel steady-state window is:
%
%       steady-state  =  [ index90_DTP  :  index90_DTR ]
%
% Adapted from DTR_10_90 and indeces10_DTRb (Bas van Niftrik / Jorn
% Fierstra, USZ), originally written for CO2/CVR stimuli.
%
% -------------------------------------------------------------------------
% INPUTS
%   BOLD4D        4D filtered fMRI data (X x Y x Z x T),  e.g. bBOLD
%   brainlag      3D voxel lag map in TR units,            e.g. lag_i10
%   DurationRA    Protocol duration vector (3 elements used):
%                   DurationRA(1) = pre-step baseline duration  (TR)
%                   DurationRA(2) = step duration               (TR)
%                   DurationRA(3) = post-step baseline duration (TR)
%   lagDiff2      Global lag offset (scalar, TR),  e.g. globalDelay_TR_i10
%   Rowscans      Scan index vector  [1, 2, ..., NrScans]
%   index90_DTP   3D map of DTP 90% index (plateau START, from computeDTP)
%   stepAmplitude 3D map of mean BOLD during step        (from computeDTP)
%
% OUTPUTS
%   DTRmap            DTR duration per voxel  (index10_DTR - index90_DTR, TR)
%   index90_DTR_map   Plateau-end index per voxel       (TR)
%   index10_DTR_map   Return-complete index per voxel   (TR)
%   Meanbaseline2map  Mean BOLD during post-step baseline per voxel
%
% -------------------------------------------------------------------------
% ASTRAN Lab - Dept of Neurosurgery, University Hospital Zurich
% Vittorio Stumpo | adapted from Bas van Niftrik / Jorn Fierstra
% -------------------------------------------------------------------------

[MtxX, MtxY, MtxZ, ~] = size(BOLD4D);
NrScans = numel(Rowscans);

% Preallocate outputs
DTRmap           = nan(MtxX, MtxY, MtxZ);
index90_DTR_map  = nan(MtxX, MtxY, MtxZ);
index10_DTR_map  = nan(MtxX, MtxY, MtxZ);
Meanbaseline2map = nan(MtxX, MtxY, MtxZ);

DurationStepchange = DurationRA(2);
Durationbaseline2  = DurationRA(3);

for ix = 1:MtxX
    for iy = 1:MtxY
        for iz = 1:MtxZ

            shift         = brainlag(ix, iy, iz);
            platStart     = index90_DTP(ix, iy, iz);    % DTP 90% = plateau START
            finalmeanstep = stepAmplitude(ix, iy, iz);  % mean BOLD during step

            % Skip invalid voxels
            if isnan(shift) || isnan(platStart) || isnan(finalmeanstep)
                continue
            end

            % Voxel-specific baseline1 duration (same convention as computeDTP)
            Durationbaseline1 = DurationRA(1) - lagDiff2 + shift;

            % Absolute TR boundaries
            stepEnd      = Durationbaseline1 + DurationStepchange;
            baseline2End = min(stepEnd + Durationbaseline2, NrScans);

            Startbaseline2 = find(Rowscans > stepEnd,      1, 'first');
            Endbaseline2   = find(Rowscans > baseline2End, 1, 'first');

            if isempty(Startbaseline2); continue; end
            if isempty(Endbaseline2)
                Endbaseline2 = NrScans;
            else
                Endbaseline2 = min(Endbaseline2, NrScans);
            end

            % Need at least 3 TRs of post-step baseline
            if (Endbaseline2 - Startbaseline2) < 3; continue; end

            fData = squeeze(BOLD4D(ix, iy, iz, :));

            try
                % --------------------------------------------------------
                % STEP 1: Initial DTR-10% index
                %   Find first TR in post-step baseline where signal has
                %   returned to within 10% of the post-step baseline level
                % --------------------------------------------------------
                BL2data = fData(Startbaseline2:Endbaseline2);
                MeanBL2 = mean(BL2data);

                BOLD10 = MeanBL2 + 0.1 * (finalmeanstep - MeanBL2);

                if BOLD10 >= MeanBL2
                    % Signal returning from above (hyperoxia-like)
                    idx10raw = find(fData(Startbaseline2:Endbaseline2) <= BOLD10, 1, 'first');
                else
                    % Signal returning from below (hypoxia-like)
                    idx10raw = find(fData(Startbaseline2:Endbaseline2) >= BOLD10, 1, 'first');
                end

                if isempty(idx10raw) || idx10raw == 0
                    IndexBOLD10dtr = Endbaseline2 - 1;
                else
                    IndexBOLD10dtr = idx10raw + Startbaseline2 - 1;
                end
                IndexBOLD10dtr = fix(min(max(IndexBOLD10dtr, Startbaseline2), Endbaseline2));

                % --------------------------------------------------------
                % STEP 2: Iterative refinement of DTR-10% (up to 15 rounds)
                %   Recomputes post-step baseline mean from current index
                %   onward, then re-finds crossing. Converges when baseline
                %   mean changes < 0.01%. Mirrors DTR_10_90 logic.
                % --------------------------------------------------------
                MeanBL2_prev = MeanBL2;

                for nIter = 1:15
                    [IndexBOLD10dtr, MeanBL2] = local_iterateDTR10(...
                        IndexBOLD10dtr, fData, finalmeanstep, Endbaseline2);

                    if abs(MeanBL2 / MeanBL2_prev - 1) < 0.0001
                        break
                    end
                    MeanBL2_prev = MeanBL2;
                end

                % --------------------------------------------------------
                % STEP 3: DTR-90% index
                %   Search backward from DTR-10% all the way to platStart
                %   (DTP90 index) to find the last TR where the signal was
                %   still at >= 90% of plateau amplitude.
                %
                %   The search extends into the step period - not just the
                %   post-step baseline - because the signal often starts
                %   returning before the stimulus ends.
                %
                %   The stimulus end plays NO role here. Decoupling DTR90
                %   from stimulus timing is the whole point. The cap
                %   min(DTR90, stimEnd) is applied later in
                %   OxygenDoubleStep_SS when computing the SS window.
                %
                %   Fallback: if no crossing found, DTR90 = stepEndIdx
                %   (signal never started returning within the step).
                % --------------------------------------------------------
                stepEndIdx = find(Rowscans > stepEnd, 1, 'first');
                if isempty(stepEndIdx)
                    stepEndIdx = fix(stepEnd);
                end
                stepEndIdx = max(fix(stepEndIdx), 1);

                BOLD90 = MeanBL2 + 0.9 * (finalmeanstep - MeanBL2);

                % Fallback: plateau extends to step end
                IndexBOLD90dtr = stepEndIdx;

                % Search backward from DTR10 to plateau start (DTP90)
                searchLimit = max(fix(platStart), 1);

                if IndexBOLD10dtr > searchLimit
                    searchVec = fData(IndexBOLD10dtr:-1:searchLimit);

                    if BOLD90 >= MeanBL2
                        idx90back = find(searchVec >= BOLD90, 1, 'first');
                    else
                        idx90back = find(searchVec <= BOLD90, 1, 'first');
                    end

                    if ~isempty(idx90back)
                        IndexBOLD90dtr = IndexBOLD10dtr - idx90back(1);
                    end
                end

                % No floor clamp at stepEndIdx - DTR90 is allowed before
                % stimulus end when the signal returns early.
                IndexBOLD90dtr = fix(IndexBOLD90dtr);
                IndexBOLD90dtr = min(IndexBOLD90dtr, IndexBOLD10dtr);

                % --------------------------------------------------------
                % STEP 4: Store results
                % --------------------------------------------------------
                DTRmap(ix, iy, iz)           = IndexBOLD10dtr - IndexBOLD90dtr;
                index90_DTR_map(ix, iy, iz)  = IndexBOLD90dtr;   % plateau END
                index10_DTR_map(ix, iy, iz)  = IndexBOLD10dtr;   % return complete
                Meanbaseline2map(ix, iy, iz) = MeanBL2;

            catch
                continue
            end

        end % iz
    end % iy

    if mod(ix, 10) == 0
        fprintf('computeDTR: slice ix = %d / %d\n', ix, MtxX);
    end

end % ix

disp('computeDTR: FINISHED');

end % main function


% =========================================================================
% LOCAL HELPER: local_iterateDTR10
%   Refines the DTR-10% index by recomputing the post-step baseline mean
%   from the current index onward. Mirrors indeces10_DTRb (Bas van Niftrik).
% =========================================================================
function [newIdx10, newMeanBL2] = local_iterateDTR10(IndexBOLD10dtr, fData, finalmeanstep, Endbaseline2)

idx10 = fix(IndexBOLD10dtr);
endBL = fix(Endbaseline2);

if idx10 >= endBL - 5
    bl2start = max(endBL - 5, 1);
else
    bl2start = idx10;
end
if bl2start >= endBL; bl2start = max(endBL - 1, 1); end

MeanBL = mean(fData(bl2start:endBL));

BOLD10 = MeanBL + 0.1 * (finalmeanstep - MeanBL);

if BOLD10 >= MeanBL
    BOLDh10 = find(fData(bl2start:endBL) <= BOLD10, 1, 'first');
else
    BOLDh10 = find(fData(bl2start:endBL) >= BOLD10, 1, 'first');
end

if isempty(BOLDh10)
    newIdx10 = endBL;
else
    newIdx10 = BOLDh10(1) + bl2start - 1;
end
newIdx10 = fix(min(max(newIdx10, bl2start), endBL));

if newIdx10 >= endBL - 5
    finalStart = max(endBL - 5, 1);
else
    finalStart = newIdx10;
end
if finalStart >= endBL; finalStart = max(endBL - 1, 1); end

newMeanBL2 = mean(fData(finalStart:endBL));

end % local_iterateDTR10
