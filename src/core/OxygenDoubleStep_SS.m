% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [HypoxiaMap_SS, HyperoxiaMap_SS, SteadyStateDurationMap_1, SteadyStateDurationMap_2, ...
          HypoxiaMap_SS_O2norm, HyperoxiaMap_SS_O2norm, DeltaPetO2_1, DeltaPetO2_2] = ...
    OxygenDoubleStep_SS(bBOLD, index90_DTP_1, index90_DTR_1, Meanbaseline_map, ...
                        index90_DTP_2, index90_DTR_2, ...
                        lag_i10, DurationRA_DTR1, DurationRA_DTR2, globalDelay_TR_i10, ...
                        UsableShiftedO2_i10)
% OxygenDoubleStep_SS  Per-voxel steady-state %%BOLD maps for double-step protocol
%
% Computes percentage BOLD change during the TRUE PLATEAU of each stimulus
% step.  The plateau window per voxel is:
%
%   START : index90_DTP
%   END   : min(index90_DTR,  stimulusEnd_per_voxel)
%
% The SS end is capped at the per-voxel stimulus end because in most
% glioblastoma voxels the BOLD signal has not yet started returning before
% the stimulus step ends. In that case index90_DTR is placed after stepEnd
% (by computeDTR's backward search), which would extend the SS window into
% the inter-step or post-step baseline.
% Only if the signal starts returning BEFORE the step ends (index90_DTR <
% stimulusEnd) is the DTR index used as the window ceiling.
%
% Additionally computes O2-normalised steady-state maps by dividing the
% raw SS %%BOLD by the achieved DeltaPetO2 for each step.  This is the
% direct O2 analogue of the CVR formula used in the lab's CO2 work:
%
%   SS_O2norm (%%BOLD/mmHg) = SS_%%BOLD / DeltaPetO2
%
% DeltaPetO2 is a per-patient scalar derived from the global
% UsableShiftedO2_i10 trace:
%
%   DeltaPetO2 = mean(O2 during baseline window)
%              - mean(O2 during step plateau window)
%
% Baseline window : scans 1 : DurationRA_DTR1(1)
% Step 1 plateau  : scans DurationRA_DTR1(1)+1 : DurationRA_DTR1(1)+DurationRA_DTR1(2)
% Step 2 plateau  : scans DurationRA_DTR2(1)+1 : DurationRA_DTR2(1)+DurationRA_DTR2(2)
% Both steps use the pre-step-1 baseline as the O2 reference, consistent
% with how Meanbaselinei10 (BOLD baseline) is defined in computeDTP.
%
% -------------------------------------------------------------------------
% INPUTS
%   bBOLD               4-D filtered fMRI data  (X x Y x Z x T)
%   index90_DTP_1       3-D map: DTP 90%% index for step 1  (plateau START, TR)
%   index90_DTR_1       3-D map: DTR 90%% index for step 1  (plateau END,   TR)
%   Meanbaseline_map    3-D map: per-voxel mean BOLD baseline
%                         (Meanbaselinei10 from computeDTP, step 1)
%   index90_DTP_2       3-D map: DTP 90%% index for step 2  (plateau START, TR)
%   index90_DTR_2       3-D map: DTR 90%% index for step 2  (plateau END,   TR)
%   lag_i10             3-D voxel lag map (TR units)
%   DurationRA_DTR1     3-element vector for step 1:
%                         [baseline1_dur, step1_dur, post_step1_dur]  (TR)
%   DurationRA_DTR2     3-element vector for step 2:
%                         [baseline2_dur, step2_dur, post_step2_dur]  (TR)
%   globalDelay_TR_i10  Global lag offset (scalar, TR)
%   UsableShiftedO2_i10 Global i10-shifted O2 trace, length = NrScans (mmHg)
%
% OUTPUTS
%   HypoxiaMap_SS             Per-voxel SS %%BOLD change, step 1
%   HyperoxiaMap_SS           Per-voxel SS %%BOLD change, step 2
%   SteadyStateDurationMap_1  Plateau duration map, step 1  (TR)
%   SteadyStateDurationMap_2  Plateau duration map, step 2  (TR)
%   HypoxiaMap_SS_O2norm      SS %%BOLD / DeltaPetO2, step 1  (%%BOLD/mmHg)
%   HyperoxiaMap_SS_O2norm    SS %%BOLD / DeltaPetO2, step 2  (%%BOLD/mmHg)
%   DeltaPetO2_1              Achieved O2 step amplitude, step 1  (mmHg)
%   DeltaPetO2_2              Achieved O2 step amplitude, step 2  (mmHg)
%
% -------------------------------------------------------------------------
% ASTRAN Lab - Dept of Neurosurgery, University Hospital Zurich
% Vittorio Stumpo / Tristan Schmidlechner
% -------------------------------------------------------------------------

[MtxX, MtxY, MtxZ, NrScans] = size(bBOLD);

HypoxiaMap_SS            = nan(MtxX, MtxY, MtxZ);
HyperoxiaMap_SS          = nan(MtxX, MtxY, MtxZ);
SteadyStateDurationMap_1 = nan(MtxX, MtxY, MtxZ);
SteadyStateDurationMap_2 = nan(MtxX, MtxY, MtxZ);

%% ------------------------------------------------------------------
%  Compute per-patient DeltaPetO2 scalars from global O2 trace
% ------------------------------------------------------------------

% Baseline: scans 1 to DurationRA_DTR1(1)
bl_end = max(1, min(fix(DurationRA_DTR1(1)), NrScans));
O2_BL  = mean(UsableShiftedO2_i10(1 : bl_end));

% Step 1 plateau window
s1_start = bl_end + 1;
s1_end   = min(bl_end + fix(DurationRA_DTR1(2)), NrScans);

% Step 2 plateau window
s2_start = fix(DurationRA_DTR2(1)) + 1;
s2_end   = min(fix(DurationRA_DTR2(1)) + fix(DurationRA_DTR2(2)), NrScans);

% Guard against empty or inverted windows
if s1_end >= s1_start
    DeltaPetO2_1 = O2_BL - mean(UsableShiftedO2_i10(s1_start : s1_end));
else
    warning('OxygenDoubleStep_SS: step 1 O2 window is empty - DeltaPetO2_1 set to NaN.');
    DeltaPetO2_1 = NaN;
end

if s2_end >= s2_start
    DeltaPetO2_2 = O2_BL - mean(UsableShiftedO2_i10(s2_start : s2_end));
else
    warning('OxygenDoubleStep_SS: step 2 O2 window is empty - DeltaPetO2_2 set to NaN.');
    DeltaPetO2_2 = NaN;
end

fprintf('OxygenDoubleStep_SS: O2_BL        = %.2f mmHg\n', O2_BL);
fprintf('OxygenDoubleStep_SS: DeltaPetO2_1 = %.2f mmHg  (scans %d:%d)\n', DeltaPetO2_1, s1_start, s1_end);
fprintf('OxygenDoubleStep_SS: DeltaPetO2_2 = %.2f mmHg  (scans %d:%d)\n', DeltaPetO2_2, s2_start, s2_end);

%% ------------------------------------------------------------------
%  Per-voxel SS %BOLD maps
% ------------------------------------------------------------------

for ix = 1:MtxX
    for iy = 1:MtxY
        for iz = 1:MtxZ

            meanBL = Meanbaseline_map(ix, iy, iz);
            shift  = lag_i10(ix, iy, iz);

            if isnan(meanBL) || meanBL == 0 || isnan(shift)
                continue
            end

            % Per-voxel stimulus end - identical formula to computeDTR:
            %   stimEnd = baseline_dur - globalDelay + voxelLag + step_dur
            stimEnd1 = DurationRA_DTR1(1) - globalDelay_TR_i10 + shift + DurationRA_DTR1(2);
            stimEnd2 = DurationRA_DTR2(1) - globalDelay_TR_i10 + shift + DurationRA_DTR2(2);

            % ------------------------------------------------------
            % Step 1
            %   START : index90_DTP_1
            %   END   : min(index90_DTR_1, stimEnd1)
            % ------------------------------------------------------
            ssStart1 = fix(index90_DTP_1(ix, iy, iz));
            ssEnd1   = min(fix(index90_DTR_1(ix, iy, iz)), fix(stimEnd1));

            if ~isnan(ssStart1) && ~isnan(ssEnd1) && ssEnd1 > ssStart1 + 1
                ssStart1 = max(ssStart1, 1);
                ssEnd1   = min(ssEnd1, NrScans);
                if ssEnd1 > ssStart1
                    meanSS1 = mean(squeeze(bBOLD(ix, iy, iz, ssStart1:ssEnd1)));
                    HypoxiaMap_SS(ix, iy, iz)            = (meanSS1 - meanBL) / meanBL * 100;
                    SteadyStateDurationMap_1(ix, iy, iz) = ssEnd1 - ssStart1;
                end
            end

            % ------------------------------------------------------
            % Step 2
            %   START : index90_DTP_2
            %   END   : min(index90_DTR_2, stimEnd2)
            % ------------------------------------------------------
            ssStart2 = fix(index90_DTP_2(ix, iy, iz));
            ssEnd2   = min(fix(index90_DTR_2(ix, iy, iz)), fix(stimEnd2));

            if ~isnan(ssStart2) && ~isnan(ssEnd2) && ssEnd2 > ssStart2 + 1
                ssStart2 = max(ssStart2, 1);
                ssEnd2   = min(ssEnd2, NrScans);
                if ssEnd2 > ssStart2
                    meanSS2 = mean(squeeze(bBOLD(ix, iy, iz, ssStart2:ssEnd2)));
                    HyperoxiaMap_SS(ix, iy, iz)          = (meanSS2 - meanBL) / meanBL * 100;
                    SteadyStateDurationMap_2(ix, iy, iz) = ssEnd2 - ssStart2;
                end
            end

        end % iz
    end % iy
end % ix

%% ------------------------------------------------------------------
%  O2-normalised SS maps  (%BOLD / mmHg)
% ------------------------------------------------------------------

if ~isnan(DeltaPetO2_1) && DeltaPetO2_1 > 0
    HypoxiaMap_SS_O2norm = HypoxiaMap_SS / DeltaPetO2_1;
else
    warning('OxygenDoubleStep_SS: DeltaPetO2_1 <= 0 or NaN - HypoxiaMap_SS_O2norm set to NaN.');
    HypoxiaMap_SS_O2norm = nan(MtxX, MtxY, MtxZ);
end

if ~isnan(DeltaPetO2_2) && DeltaPetO2_2 > 0
    HyperoxiaMap_SS_O2norm = HyperoxiaMap_SS / DeltaPetO2_2;
else
    warning('OxygenDoubleStep_SS: DeltaPetO2_2 <= 0 or NaN - HyperoxiaMap_SS_O2norm set to NaN.');
    HyperoxiaMap_SS_O2norm = nan(MtxX, MtxY, MtxZ);
end

disp('OxygenDoubleStep_SS: FINISHED');

end
